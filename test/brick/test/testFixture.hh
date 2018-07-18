/**
***************************************************************************
* @file brick/test/testFixture.hh
*
* Header file declaring TestFixture class.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_TEST_TESTFIXTURE_HH
#define BRICK_TEST_TESTFIXTURE_HH

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <brick/test/runnableObject.hh>
#include <brick/test/testException.hh>
#include <brick/test/testMacros.hh>

namespace brick {

  namespace test {

    /**
     ** The TestFixture class helps with unit testing by coordinating
     ** the execution of a test suite.  To use this class, you should
     ** create a derived class using the following slightly scary
     ** inheritance:
     **
     **   class MyTestClass : public brick::TestFixture<MyTestClass>;
     **
     ** The derived class should define individual tests as member
     ** functions:
     **
     **   void MyTestClass::testFunction0() {...};
     **
     ** These tests should be registered with the TestFixture, probably
     ** in the constructor of the derived class.  You can do this the
     ** easy way, by using the provided macro:
     **
     **   BRICK_TEST_REGISTER_MEMBER(testFunction0);
     **
     ** Note that this macro will fail if you call it from any context
     ** that is not a constructor or member function of the derived
     ** class.
     **
     ** Or you can do it the hard way by calling the registerTest member
     ** function directly:
     **
     **   this->registerTest("testFunction0", &MyTestClass::testFunction0);
     **
     ** Within the tests, use the BRICK_TEST_ASSERT macros on pass/fail
     ** conditions:
     **
     **   BRICK_TEST_ASSERT(residual < m_errorThreshold);
     **
     **   BRICK_TEST_ASSERT_EXCEPTION(ValueException, functionWhichThrows());
     **
     ** Then you can run the tests by calling the run() method, which
     ** your derived class will have inherited from TestFixture<>.
     **/
    template <class FixtureType>
    class TestFixture
      : public RunnableObject
    {
    public:

      /* ============== Public Typedefs ============== */

      /**
       * The TestFixtureType typedef is used by helper macros, and may
       * eventually be useful for other things, too.
       */
      typedef FixtureType TestFixtureType;

      /**
       * The TestFunctionPtr typedef is used to conveniently work with
       * member functions of the derived class, which will be run as
       * individual tests.
       */
      typedef void (FixtureType::* TestFunctionPtr)();


      /* ============== Public Member Functions ============== */

      /**
       * This constructor sets the name of the group of tests run by
       * this test fixture.
       *
       * @param testFixtureName This argument specifies the name by
       * which this test fixture will be referred to in status output.
       */
      explicit
      TestFixture(const std::string& testFixtureName);


      /**
       * The destructor cleans up any resources and destroys the class
       * instance.
       */
      virtual
      ~TestFixture();


      /**
       * This member function is used to indicate which member functions
       * of the subclass should be run as tests.  Call it once for each
       * test.
       *
       * @param testName This argument specifies the name of the test.
       * It will be used to identify the test in text output, and it
       * will be passed to both the setUp() member function (before the
       * test is run) and the tearDown() member function (after the test
       * is run).
       *
       * @param testFunctionPtr This argument should be a pointer to the
       * member function which should be run as a test.
       */
      virtual void
      registerTest(const std::string& testName,
                   TestFunctionPtr testFunctionPtr);

      /**
       * This member function runs all of the registered tests, keeps
       * track of the results, and prints any diagnostic output, finally
       * returning a bool indicating success (true) or failure (false).
       *
       * @return The return value will be true if all of the registered
       * tests pass, false otherwise.
       */
      virtual bool
      run();

    protected:

      /* ============== Protected Member Functions ============== */

      /**
       * This protected member function prints a message indicating that
       * the test fixture has run all of its tests.  If a subclass
       * overrides TestFixture::run(), then this method should generally
       * be called at the end of the run() method.
       */
      virtual void
      announceTestFinish();

      /**
       * This protected member function prints a message indicating that
       * the test fixture is about to run its tests.  If a subclass
       * overrides TestFixture::run(), then this method should generally
       * be called at the beginning of the run() method.
       */
      virtual void
      announceTestStart();


      /**
       * This protected member function builds a diagnostic string
       * describing a test failure.  Derived classes can override this
       * to customize test failure output.
       *
       * @param failureIndex This argument should be set to the number of
       * the failure.  It will be included in the test message.  Increment
       * this number each time you call buildFailureMessage().
       *
       * @param testName This argument specifies the name of the failed
       * test.
       *
       * @param failureType This argument specifies the type of
       * failure.  Suggestions are "Failed test" for tests which fail,
       * and "Error in test" for tests which crash.
       *
       * @param whatMessage This argument is set to the "what()"
       * output of the caught exception, and should usually provide
       * information about the nature of the failure.
       *
       * @return The return value is a formatted string describing the
       * error.  The default format looks something like this:
       * "1) Failed test foo::bar ..."
       */
      virtual std::string
      buildFailureMessage(size_t failureIndex,
                          const std::string& testName,
                          const std::string& failureType,
                          const std::string& whatMessage);


      /**
       * This protected member function prints the final output
       * indicating how many tests passed/failed, etc.
       *
       * @param errorCount This argument specifies how many tests
       * aborted by throwing exceptions.
       *
       * @param testMessages This argument is a vector of failure
       * messages, one for each test that didn't pass.
       */
      void
      printTestStatistics(int errorCount,
                          const std::vector<std::string>& testMessages);

      /**
       * This protected member function is called immediately before
       * each test is run.  The argument will be set to the name of the
       * test which is about to be run.
       *
       * @param testName The argument will be set to the name of the
       * test which is about to be run, allowing the setUp() method of
       * the derived class to modify its behavior for each test.
       */
      virtual void
      setUp(const std::string&) {}

      /**
       * This protected member function is called immediately after
       * each test is run.  The argument will be set to the name of the
       * test which was just run.
       *
       * @param testName The argument will be set to the name of the
       * test which was just run, allowing the tearDown() method of
       * the derived class to modify its behavior for each test.
       */
      virtual void
      tearDown(const std::string&) {}


      /* ============== Protected Member Variables ============== */

      /**
       * This variable controls at what character position text output
       * will be wrapped.
       */
      int m_textOutputLineLength;

      /**
       * This variable controls the amount of text output produced while
       * running tests.  Lower numbers mean less output.
       */
      int m_verbosity;

    private:

      /* ============== Private Member Variables ============== */


      /* ============== Private Member Variables ============== */

      std::string m_testFixtureName;
      std::vector<TestFunctionPtr> m_testFunctionPtrVector;
      std::vector<std::string> m_testNameVector;
    };

  } // namespace test

} // namespace brick


/* =================================================================
 * Member function definitions follow.
 * This would be a .cc file if TestFixture weren't a template.
 * ================================================================= */

namespace brick {

  namespace test {

    // This constructor sets the name of the group of tests run by
    // this test fixture.
    template <class FixtureType>
    TestFixture<FixtureType>::
    TestFixture(const std::string& testFixtureName)
      : RunnableObject(),
	m_textOutputLineLength(75),
        m_verbosity(1),
        m_testFixtureName(testFixtureName),
        m_testFunctionPtrVector(),
        m_testNameVector()
    {
      // Empty
    }


    // The destructor cleans up any resources and destroys the class
    // instance.
    template <class FixtureType>
    TestFixture<FixtureType>::
    ~TestFixture()
    {
      // Empty
    }


    // This member function is used to indicate which member functions
    // of the subclass should be run as tests.
    template <class FixtureType>
    void
    TestFixture<FixtureType>::
    registerTest(const std::string& testName,
                 typename TestFixture<FixtureType>::TestFunctionPtr
                 testFunctionPtr)
    {
      m_testNameVector.push_back(testName);
      m_testFunctionPtrVector.push_back(testFunctionPtr);
    }


    // This member function runs all of the registered tests, keeps
    // track of the results, and prints any diagnostic output, finally
    // returning a bool indicating success (true) or failure (false).
    template <class FixtureType>
    bool
    TestFixture<FixtureType>::
    run()
    {
      // Print the starting banner.
      this->announceTestStart();

      // This vector will hold the a message for each failed test.
      std::vector<std::string> testMessages;

      // This counter will reflect how many tests fail by throwing
      // non-test exceptions.
      int errorCount = 0;

      // Now run each test in turn.
      for(size_t index = 0; index < m_testFunctionPtrVector.size(); ++index) {
        // Run the setUp() member function to initialize any test data.
        try {
          this->setUp(m_testNameVector[index]);
        } catch(const std::exception&) {
          // Is there anything useful to do here other than just re-throw
          // the exception?
          throw;
        }

        // Extract the function pointer.
        TestFunctionPtr testFunctionPtr = m_testFunctionPtrVector[index];
        try {
          // Run the test.
          FixtureType* subclassThis = dynamic_cast<FixtureType*>(this);
          (subclassThis->*testFunctionPtr)();
          // Successful completion.  Print a happy ".".
          std::cout << "." << std::flush;

        } catch(const TestException& caughtException) {
          // Failure!  The test threw a TestException, the purpose of which
          // is to indicate a failed test.
          // Write status output.
          std::cout << "F" << std::flush;

          // And save information about the failure.
          std::string message = this->buildFailureMessage(
            testMessages.size() + 1, m_testNameVector[index], "Failed test",
            caughtException.what());
          testMessages.push_back(message);

        } catch(const std::exception& caughtException) {
          // Error!  The test threw an std::exception besides TestException.
          // Write status output.
          std::cout << "E" << std::flush;

          // Increment the count of errors.
          ++errorCount;

          // And save information about the error.
          std::string message = this->buildFailureMessage(
            testMessages.size() + 1, m_testNameVector[index], "Error in test",
            caughtException.what());
          testMessages.push_back(message);

        } catch(...) {
          // Error!  The test threw something, but who knows what!
          // Write status output.
          std::cout << "E" << std::flush;

          // Increment the count of errors.
          ++errorCount;

          // And save information about the error.
          std::string message = this->buildFailureMessage(
            testMessages.size() + 1, m_testNameVector[index], "Error in test",
            "Unrecognized exception");
          testMessages.push_back(message);
        }

        // Run the tearDown() member function to clean up.
        try {
          this->tearDown(m_testNameVector[index]);
        } catch(const std::exception&) {
          // Is there anything useful to do here other than just re-throw
          // the exception?
          throw;
        }

        // We respect 80 column terminals by printing a carriage return
        // every once in a while.
        if((index % m_textOutputLineLength == 0)
           && (index != 0)) {
          std::cout << std::endl;
        }
      }

      // Print a summary of what happened.
      this->printTestStatistics(errorCount, testMessages);

      // Print the ending banner and return.
      this->announceTestFinish();
      return (testMessages.size() == 0);
    }


    // This protected member function prints a message indicating that
    // the test fixture has run all of its tests.
    template <class FixtureType>
    void
    TestFixture<FixtureType>::
    announceTestFinish()
    {
      if(m_verbosity >= 3) {
        std::cout << "/////////////////////////////////////////////////////\n"
                  << "// Completed test: " << m_testFixtureName << "\n"
                  << "/////////////////////////////////////////////////////\n"
                  << std::endl;
      }
    }


    // This protected member function prints a message indicating that
    // the test fixture is about to run its tests.
    template <class FixtureType>
    void
    TestFixture<FixtureType>::
    announceTestStart()
    {
      if(m_verbosity >= 1) {
        std::cout << "/////////////////////////////////////////////////////\n"
                  << "// Starting test: " << m_testFixtureName << "\n"
                  << "/////////////////////////////////////////////////////\n"
                  << std::endl;
      }
    }


    // This protected member function builds a diagnostic string
    // describing a test failure.
    template <class FixtureType>
    std::string
    TestFixture<FixtureType>::
    buildFailureMessage(size_t failureIndex,
                        const std::string& testName,
                        const std::string& failureType,
                        const std::string& whatMessage)
    {
      // We'll build the message using stream IO.
      std::ostringstream messageBuffer;

      // First format up the "1) " part.
      messageBuffer << failureIndex << ") ";

      // We'll want to know how many spaces the string takes up so far
      // so we can format better later.
      size_t indentSize = messageBuffer.str().size();

      // Add the failure description and test name.
      messageBuffer << failureType << " " << m_testFixtureName << "::"
                    << testName;

      // Now add the exception output, blindly breaking the text to fit
      // the (usually 80 column) screen.
      std::string::size_type currentIndex = 0;
      std::string::size_type chunkSize = m_textOutputLineLength - indentSize;

      std::string whitespace = " \t\n";
      while(currentIndex < whatMessage.size()) {
        // Carriage return and indent as appropriate.
        messageBuffer << "\n";
        for(size_t index = 0; index < indentSize; ++index) {
          messageBuffer << " ";
        }

        // If the line is too long to fit without a line break, then
        // search for a good whitespace position at which to break the
        // line.  When this loop terminates, whitespaceIndex should
        // point to the whitespace that comes closest to the ideal line
        // break position (chunkSize) without exceeding it.  If there's
        // no whitespace within the target line length, then
        // whitespaceIndex will be equal to currentIndex.  In either
        // case, nextWhitespaceIndex will point to the nearest
        // whitespace which exceeds the ideal line break length, or npos
        // if there's no whitespace at all.
        std::string::size_type whitespaceIndex = currentIndex;
        std::string::size_type nextWhitespaceIndex = currentIndex;
        if(whatMessage.size() - currentIndex > chunkSize) {
          while(nextWhitespaceIndex <= currentIndex + chunkSize) {
            whitespaceIndex = nextWhitespaceIndex;
            nextWhitespaceIndex =
              whatMessage.find_first_not_of(whitespace, whitespaceIndex);
            if(nextWhitespaceIndex == std::string::npos) {
              break;
            }
            nextWhitespaceIndex =
              whatMessage.find_first_of(whitespace, nextWhitespaceIndex);
          }
        }

        // Write the next chunk of the string and move to the next chunk.
        if(whitespaceIndex != currentIndex) {
          std::string::size_type localChunkSize =
            whitespaceIndex - currentIndex;
          messageBuffer << whatMessage.substr(currentIndex, localChunkSize);
          currentIndex += localChunkSize + 1;
        } else {
          messageBuffer << whatMessage.substr(currentIndex, chunkSize);
          currentIndex += chunkSize;
        }
      }

      // Done with formatting.
      return messageBuffer.str();
    }


    // This protected member function prints the final output
    // indicating how many tests passed/failed, etc.
    template <class FixtureType>
    void
    TestFixture<FixtureType>::
    printTestStatistics(int errorCount,
                        const std::vector<std::string>& testMessages)
    {
      // Print statistics from the test run.
      std::cout << "\n\n";
      std::cout << "Test Results:\n"
                << "Run: " << m_testFunctionPtrVector.size()
                << "  \tFailures: " << testMessages.size() - errorCount
                << "  \tErrors: " << errorCount << "\n" << std::endl;

      // Print accumulated error messages.
      for(size_t index = 0; index < testMessages.size(); ++index) {
        std::cout << testMessages[index] << std::endl;
        std::cout << std::endl;
      }
    }

  } // namespace test

} // namespace brick

#endif /* #ifdef BRICK_TEST_TESTFIXTURE_HH */
