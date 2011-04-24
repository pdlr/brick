/**
***************************************************************************
* @file brick/utilties/test/lockFileTest.cc
* 
* Source file defining tests for lock file handling.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <fstream>
#include <sstream>
#include <string>

#include <brick/portability/timeUtilities.hh>
#include <brick/test/testFixture.hh>
#include <brick/utilities/lockFile.hh>
#include <brick/utilities/path.hh>

using brick::test::TestFixture;

namespace brick {

  namespace utilities {
    
    class LockFileTest : public TestFixture<LockFileTest> {

    public:

      LockFileTest();
      ~LockFileTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Tests of member functions.
      void testConstructor__string();
      void testConstructor__string__string();
      void testDestructor();
      void testIsValid();

    private:

      std::string
      getInvalidLockFileName();

    
      std::string
      getLockFileName();
    
    }; // class LockFileTest


    /* ============== Member Function Definititions ============== */

    LockFileTest::
    LockFileTest()
      : TestFixture<LockFileTest>("LockFileTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testConstructor__string);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__string__string);
      BRICK_TEST_REGISTER_MEMBER(testDestructor);
      BRICK_TEST_REGISTER_MEMBER(testIsValid);
    }


    void
    LockFileTest::
    testConstructor__string()
    {
      // Create a lock file names.
      std::string lockFileName = this->getLockFileName();
      std::string invalidLockFileName = this->getInvalidLockFileName();

      if(isExistingPath(lockFileName) || isExistingPath(invalidLockFileName)) {
        BRICK_THROW(brick::common::IOException,
                    "LockFileTest::testConstructor__string()",
                    "Lock files already exist.");
      }

      {
        LockFile lockFile0(lockFileName);
        LockFile lockFile1(lockFileName);
        LockFile lockFile2(lockFileName);
        LockFile lockFile3(invalidLockFileName);
        LockFile lockFile4(invalidLockFileName);

        BRICK_TEST_ASSERT(isRegularFile(lockFileName));
        BRICK_TEST_ASSERT(!isRegularFile(invalidLockFileName));

        BRICK_TEST_ASSERT(lockFile0.isValid());
        BRICK_TEST_ASSERT(!(lockFile1.isValid()));
        BRICK_TEST_ASSERT(!(lockFile2.isValid()));
        BRICK_TEST_ASSERT(!(lockFile3.isValid()));
        BRICK_TEST_ASSERT(!(lockFile4.isValid()));
      }

      BRICK_TEST_ASSERT(!isRegularFile(lockFileName));

      {
        LockFile lockFile0(lockFileName);
        LockFile lockFile1(lockFileName);
        BRICK_TEST_ASSERT(lockFile0.isValid());
        BRICK_TEST_ASSERT(!(lockFile1.isValid()));
      }

      BRICK_TEST_ASSERT(!isRegularFile(lockFileName));
    }

  
    void
    LockFileTest::
    testConstructor__string__string()
    {
      // Create a lock file names.
      std::string lockFileName = this->getLockFileName();
      std::string invalidLockFileName = this->getInvalidLockFileName();

      // Crate string for lock file contents.
      std::string contents = "It's not unusual\nto be loved by anyone.";
    
      if(isExistingPath(lockFileName) || isExistingPath(invalidLockFileName)) {
        BRICK_THROW(brick::common::IOException,
                    "LockFileTest::testConstructor__string()",
                    "Lock files already exist.");
      }

      {
        LockFile lockFile0(lockFileName, contents);
        LockFile lockFile1(lockFileName, contents);
        LockFile lockFile2(invalidLockFileName, contents);

        BRICK_TEST_ASSERT(isRegularFile(lockFileName));
        BRICK_TEST_ASSERT(!isRegularFile(invalidLockFileName));

        BRICK_TEST_ASSERT(lockFile0.isValid());
        BRICK_TEST_ASSERT(!(lockFile1.isValid()));
        BRICK_TEST_ASSERT(!(lockFile2.isValid()));

        std::ifstream inputFileStream(lockFileName.c_str());
        std::istringstream inputStringStream(contents);
        while(inputStringStream.good()) {
          BRICK_TEST_ASSERT(inputFileStream.good());
          BRICK_TEST_ASSERT(inputFileStream.get() == inputStringStream.get());
        }
        BRICK_TEST_ASSERT(!inputFileStream.good());
      }

      BRICK_TEST_ASSERT(!isRegularFile(lockFileName));
    }


    void
    LockFileTest::
    testDestructor()
    {
      // Tested by constructor tests.
    }

  
    void
    LockFileTest::
    testIsValid()
    {
      // Tested by constructor tests.
    }

  
    std::string
    LockFileTest::
    getInvalidLockFileName()
    {
      std::string directoryName;
      int index0 = 0;
      while(1) {
        // We don't use portable file naming conventions here, but it
        // doesn't matter because we're trying to create a nonexistant
        // directory name anyway.
        std::ostringstream directoryNameStream;
        directoryNameStream << "/var/tmp/nonExistant" << index0;
        directoryName = directoryNameStream.str();
        if(!isDirectory(directoryName)) {
          break;
        }
        ++index0;
      }
      // We _do_ use portable file naming here to avoid accidentally
      // creating a valid filename.
      std::string fileName = joinPath(directoryName, "foo.tmp");
      return fileName;
    }

  
    std::string
    LockFileTest::
    getLockFileName()
    {
      std::string lockFileName;
      while(1) {
        std::ostringstream lockFileNameStream;
        lockFileNameStream << "lockFile" << portability::getCurrentTime()
                           << ".tmp";
        lockFileName = lockFileNameStream.str();
        if(!isExistingPath(lockFileName)) {
          break;
        }
      }
      return lockFileName;
    }
  
  } // namespace utilities
  
} // namespace brick


#ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION

int main(int argc, char** argv)
{
  brick::utilities::LockFileTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */

namespace {

  brick::utilities::LockFileTest currentTest;
  
}

#endif /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */
