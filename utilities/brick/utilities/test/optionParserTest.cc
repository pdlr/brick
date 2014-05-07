/**
***************************************************************************
* @file brick/utilities/test/optionParserTest.cc
* 
* Source file defining tests command line parsing.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <string>

#include <brick/common/mathFunctions.hh>
#include <brick/utilities/optionParser.hh>
#include <brick/test/testFixture.hh>

using brick::test::TestFixture;

namespace brick {

  namespace utilities {

    class OptionParserTest : public TestFixture<OptionParserTest> {

    public:

      OptionParserTest();
      ~OptionParserTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      void testConstructor();
      void testDestructor();
      void testAddOption();
      void testAddOptionWithValue();
      void testAddPositionalArgument();
      void testGetCount();
      void testGetExtraPositionalArguments();
      void testGetOptionsDescription();
      void testGetUsage();
      void testGetValue();
      void testParseCommandLine();

    private:

    }; // class OptionParserTest


    /* ============== Member Function Definititions ============== */

    OptionParserTest::
    OptionParserTest()
      : TestFixture<OptionParserTest>("OptionParserTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testConstructor);
      BRICK_TEST_REGISTER_MEMBER(testDestructor);
      BRICK_TEST_REGISTER_MEMBER(testAddOption);
      BRICK_TEST_REGISTER_MEMBER(testAddOptionWithValue);
      BRICK_TEST_REGISTER_MEMBER(testAddPositionalArgument);
      BRICK_TEST_REGISTER_MEMBER(testGetCount);
      BRICK_TEST_REGISTER_MEMBER(testGetExtraPositionalArguments);
      BRICK_TEST_REGISTER_MEMBER(testGetOptionsDescription);
      BRICK_TEST_REGISTER_MEMBER(testGetUsage);
      BRICK_TEST_REGISTER_MEMBER(testGetValue);
      BRICK_TEST_REGISTER_MEMBER(testParseCommandLine);
    }


    void
    OptionParserTest::
    testConstructor()
    {
      int argc0 = 6;
      const char* argv0[] = {"foo", "-b", "10", "--bar=Ten", "--baz", "baz"};
      int argc1 = 2;
      const char* argv1[] = {"foo", "bar"};
      int argc2 = 3;
      const char* argv2[] = {"foo", "-blr", "-prl"};

      {
        // Default parser.
        OptionParser optionParser0;

        // Make sure pre-parsing behavior is sane.
        BRICK_TEST_ASSERT_EXCEPTION(brick::common::ValueException,
                                    optionParser0.getCount("Foo"));
        BRICK_TEST_ASSERT_EXCEPTION(brick::common::ValueException,
                                    optionParser0.getValue("Foo"));

        // This should throw because of unexpected options.
        BRICK_TEST_ASSERT_EXCEPTION(
          brick::common::IOException,
          optionParser0.parseCommandLine(argc0, argv0));
      }

    
      {
        // Test "allowExtraArguments".
        OptionParser optionParser0;
        OptionParser optionParser1(false);
        optionParser0.parseCommandLine(argc1, argv1);
        std::vector<std::string> extraArguments = 
          optionParser0.getExtraPositionalArguments();
        BRICK_TEST_ASSERT(extraArguments.size() == 1);
        BRICK_TEST_ASSERT(extraArguments[0] == "bar");
        BRICK_TEST_ASSERT_EXCEPTION(
          brick::common::IOException,
          optionParser1.parseCommandLine(argc1, argv1));
      }

    
      {
        // Test allowStackedShortOptions.

        // This parser should allow stacked arguments in which, for
        // example, "-blr" is equivalent to "-b -l -r".
        OptionParser optionParser0(true);

        // This parser should not.
        OptionParser optionParser1(true, false);

        optionParser0.addOption("Bar", "-b", "--bar", "Specify bar.");
        optionParser0.addOption("Lar", "-l", "--lar", "Specify lar.");
        optionParser0.addOption("Rar", "-r", "--rar", "Specify rar.");
        optionParser0.addOptionWithValue(
          "Par", "-p", "--par", "", "Specify par.");
        optionParser1.addOption("Bar", "-b", "--bar", "Specify bar.");
        optionParser1.addOption("Lar", "-l", "--lar", "Specify lar.");
        optionParser1.addOption("Rar", "-r", "--rar", "Specify rar.");
        optionParser1.addOptionWithValue(
          "Par", "-p", "--par", "", "Specify par.");
      
        optionParser0.parseCommandLine(argc2, argv2);
        BRICK_TEST_ASSERT(optionParser0.getCount("Bar") == 1);
        BRICK_TEST_ASSERT(optionParser0.getCount("Lar") == 1);
        BRICK_TEST_ASSERT(optionParser0.getCount("Rar") == 1);
        BRICK_TEST_ASSERT(optionParser0.getCount("Par") == 1);
        BRICK_TEST_ASSERT(optionParser0.getValue("Par") == "rl");
      
        BRICK_TEST_ASSERT_EXCEPTION(
          brick::common::IOException,
          optionParser1.parseCommandLine(argc2, argv2));
      }

    
      {
        // Test allowOptionishArguments.
        // This parser should allow positional arguments that start with "-".
        OptionParser optionParser1(true, true, true);

        // We expect the same pre-parsing behavior.
        BRICK_TEST_ASSERT_EXCEPTION(brick::common::ValueException,
                                    optionParser1.getCount("Foo"));
        BRICK_TEST_ASSERT_EXCEPTION(brick::common::ValueException,
                                    optionParser1.getValue("Foo"));

        // All arguments should be interpreted as positional.
        optionParser1.parseCommandLine(argc0, argv0);
        std::vector<std::string> extraArguments =
          optionParser1.getExtraPositionalArguments();
        BRICK_TEST_ASSERT(int(extraArguments.size()) == argc0 - 1);
        for(size_t argIndex = 0; argIndex < extraArguments.size(); ++argIndex) {
          BRICK_TEST_ASSERT(extraArguments[argIndex] == argv0[argIndex + 1]);
        }
      }
    }

  
    void
    OptionParserTest::
    testDestructor()
    {
      // No explicit test for destructor.
    }

  
    void
    OptionParserTest::
    testAddOption()
    {
      int argc0 = 6;
      const char* argv0[] = {"foo", "-b", "10", "--bar=Ten", "--baz", "baz"};
      int argc1 = 6;
      const char* argv1[] = {"foo", "-b", "10", "--bar", "--baz", "baz"};
      int argc2 = 6;
      const char* argv2[] = {"foo", "-bbl", "10", "--bar", "--baz", "baz"};
    
      OptionParser optionParser0;
      optionParser0.addOption("BAR", "-b", "--bar", "Set the bar.");
      optionParser0.addOption("BAZ", "", "--baz", "Set the baz.");

      // Should throw because of conflict with previously set BAZ
      // option.
      BRICK_TEST_ASSERT_EXCEPTION(
        brick::common::StateException,
        optionParser0.addOption("BAZ", "", "--baz", "Set the baz."));

      // Should throw because --bar takes no argument.
      BRICK_TEST_ASSERT_EXCEPTION(
        brick::common::IOException,
        optionParser0.parseCommandLine(argc0, argv0));

      // Should succeed.
      optionParser0.parseCommandLine(argc1, argv1);
      BRICK_TEST_ASSERT(optionParser0.getCount("BAR") == 2);
      BRICK_TEST_ASSERT(optionParser0.getCount("BAZ") == 1);

      // Test parsing short-style options which are followed immediately
      // by values, like "-p10".
      OptionParser optionParser1(true);
      optionParser1.addOption("BAR", "-b", "--bar", "Set the bar.");
      optionParser1.addOption("BAZ", "", "--baz", "Set the baz.");

      // Should throw because -bbl parses to "-b bl" and -b does not
      // take an argument.
      BRICK_TEST_ASSERT_EXCEPTION(
        brick::common::IOException,
        optionParser1.parseCommandLine(argc2, argv2));
    }

  
    void
    OptionParserTest::
    testAddOptionWithValue()
    {
      int argc0 = 6;
      const char* argv0[] = {"foo", "-b", "10", "--bar=Ten", "--baz", "baz"};
      int argc1 = 6;
      const char* argv1[] = {"foo", "-bbl", "10", "--bar=Ten", "--baz", "baz"};
      int argc2 = 1;
      const char* argv2[] = {"foo"};
    
      OptionParser optionParser0;
      optionParser0.addOptionWithValue(
        "BAR", "-b", "--bar", "7", "Set the bar.");
      optionParser0.addOptionWithValue(
        "BAZ", "", "--baz", "6", "Set the baz.");
      optionParser0.addOptionWithValue(
        "PRY", "-p", "--pry", "--", "Set the pry.");

      optionParser0.parseCommandLine(argc0, argv0);
      BRICK_TEST_ASSERT(optionParser0.getCount("BAR") == 2);
      BRICK_TEST_ASSERT(optionParser0.getCount("BAZ") == 1);
      BRICK_TEST_ASSERT(optionParser0.getCount("PRY") == 0);
      BRICK_TEST_ASSERT(optionParser0.getValue("BAR") == "Ten");
      BRICK_TEST_ASSERT(optionParser0.getValue("BAZ") == "baz");
      BRICK_TEST_ASSERT(optionParser0.getValue("PRY") == "--");

      // Check that multiple occurrences of the same option are saved.
      BRICK_TEST_ASSERT(optionParser0.getValue("BAR", 0) == "10");
      BRICK_TEST_ASSERT(optionParser0.getValue("BAR", 1) == "Ten");
    
      OptionParser optionParser1;
      optionParser1.addOptionWithValue(
        "BAR", "", "--bar", "", "Set the bar.");
      optionParser1.addOptionWithValue(
        "BAZ", "", "--baz", "6", "Set the baz.");
      optionParser1.addOptionWithValue(
        "BASTE", "-b", "--baste", "7", "Set the baste.");
      optionParser1.parseCommandLine(argc0, argv0);
 
      BRICK_TEST_ASSERT(optionParser1.getCount("BAR") == 1);
      BRICK_TEST_ASSERT(optionParser1.getCount("BAZ") == 1);
      BRICK_TEST_ASSERT(optionParser1.getCount("BASTE") == 1);
      BRICK_TEST_ASSERT(optionParser1.getValue("BAR") == "Ten");
      BRICK_TEST_ASSERT(optionParser1.getValue("BAZ") == "baz");
      BRICK_TEST_ASSERT(optionParser1.getValue("BASTE") == "10");


      // Test parsing short-style options which are followed immediately
      // by values, like "-p10".
      optionParser1.parseCommandLine(argc1, argv1);
      std::vector<std::string> extraArguments =
        optionParser1.getExtraPositionalArguments();
 
      BRICK_TEST_ASSERT(optionParser1.getCount("BAR") == 1);
      BRICK_TEST_ASSERT(optionParser1.getCount("BAZ") == 1);
      BRICK_TEST_ASSERT(optionParser1.getCount("BASTE") == 1);
      BRICK_TEST_ASSERT(optionParser1.getValue("BAR") == "Ten");
      BRICK_TEST_ASSERT(optionParser1.getValue("BAZ") == "baz");
      BRICK_TEST_ASSERT(optionParser1.getValue("BASTE") == "bl");
      BRICK_TEST_ASSERT(extraArguments.size() == 1);
      BRICK_TEST_ASSERT(extraArguments[0] == "10");

      // Check that non-string default arguments are handled OK.
      OptionParser optionParser2;
      optionParser2.addOptionWithValue(
        "BAR", "-b", "--bar", 7, "Set the bar.");
      optionParser2.addOptionWithValue(
        "BAZ", "", "--baz", 6.5, "Set the baz.");
      optionParser2.addOptionWithValue(
        "PRY", "-p", "--pry", "--", "Set the pry.");

      optionParser0.parseCommandLine(argc2, argv2);
      BRICK_TEST_ASSERT(optionParser2.convertValue<int>("BAR") == 7);
      BRICK_TEST_ASSERT(
        brick::common::absoluteValue(optionParser2.convertValue<double>("BAZ")
                                     - 6.5) < 1.0e-9);
    }

  
    void
    OptionParserTest::
    testAddPositionalArgument()
    {
      int argc0 = 6;
      const char* argv0[] = {"foo", "-b", "10", "--bar", "--baz", "baz"};
    
      // Test parsing with required and optional arguments.
      OptionParser optionParser0;
      optionParser0.addOption("BAR", "-b", "--bar", "Set the bar.");
      optionParser0.addOption("BAZ", "", "--baz", "Set the baz.");
      optionParser0.addPositionalArgument("FirstPos", "The first pos", true);
      optionParser0.addPositionalArgument("SecondPos", "The second pos", true);
      optionParser0.addPositionalArgument("ThirdPos", "The third pos");
      optionParser0.addPositionalArgument(
        "FourthPos", "The fourth pos", false, "dummyValue");

      optionParser0.parseCommandLine(argc0, argv0);
      BRICK_TEST_ASSERT(optionParser0.getValue("FirstPos") == "10");
      BRICK_TEST_ASSERT(optionParser0.getValue("SecondPos") == "baz");
      BRICK_TEST_ASSERT(optionParser0.getValue("ThirdPos") == "");
      BRICK_TEST_ASSERT(optionParser0.getValue("FourthPos") == "dummyValue");

      // Test parsing with required positional arguments.
      OptionParser optionParser1;
      optionParser1.addOption("BAR", "-b", "--bar", "Set the bar.");
      optionParser1.addOption("BAZ", "", "--baz", "Set the baz.");
      optionParser1.addPositionalArgument("FirstPos", "The first pos", true);
      optionParser1.addPositionalArgument("SecondPos", "The second pos", true);
      optionParser1.addPositionalArgument("ThirdPos", "The third pos", true);

      BRICK_TEST_ASSERT_EXCEPTION(
        brick::common::IOException,
        optionParser1.parseCommandLine(argc0, argv0));
    }

  
    void
    OptionParserTest::
    testGetCount()
    {
      // Already tested in testAddOption().
    }

  
    void
    OptionParserTest::
    testGetExtraPositionalArguments()
    {
      // Already tested in testConstructor().
    }

  
    void
    OptionParserTest::
    testGetOptionsDescription()
    {
      // No explicit test.
    }

  
    void
    OptionParserTest::
    testGetUsage()
    {
      // No explicit test.
    }

  
    void
    OptionParserTest::
    testGetValue()
    {
      // Already tested in testAddOptionWithValue().
    }

  
    void
    OptionParserTest::
    testParseCommandLine()
    {
      // Already tested in testAddOptionWithValue().
    }

  } // namespace utilities
  
} // namespace brick


#ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION

int main(int argc, char** argv)
{
  brick::utilities::OptionParserTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */

namespace {

  brick::utilities::OptionParserTest currentTest;

}

#endif /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */
