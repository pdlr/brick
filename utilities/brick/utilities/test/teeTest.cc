/**
***************************************************************************
* @file brick/utilities/test/teeTest.cc
*
* Source file defining tests for the Tee class.
*
* Copyright (C) 2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <fstream>
#include <string>

#include <brick/utilities/tee.hh>
#include <brick/test/testFixture.hh>

using brick::test::TestFixture;

namespace brick {

  namespace utilities {

    class TeeTest : public TestFixture<TeeTest> {

    public:

      TeeTest();
      ~TeeTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Tests of member functions.
      void testTee();

    }; // class TeeTest


    /* ============== Member Function Definititions ============== */

    TeeTest::
    TeeTest()
      : TestFixture<TeeTest>("TeeTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testTee);
    }


    void
    TeeTest::
    testTee()
    {
#if 0
      TeeBuffer<char> teeBuffer;
      teeBuffer.add(std::cout.rdbuf());
      teeBuffer.add(std::cerr.rdbuf());
      std::ostream testStream(&teeBuffer);
      testStream << "This is a test" << std::endl;
#else
      Tee<char> tee;
      tee.add(std::cout);
      tee.add(std::cerr);
      tee << "This is a test" << std::endl;
#endif
    }

  } // namespace utilities
  
} // namespace brick

// xxx
#ifndef BRICK_TEST_NO_AUTOMATIC_REGISTRATION
#define BRICK_TEST_NO_AUTOMATIC_REGISTRATION 1
#endif

#ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION

int main(int /* argc */, char** /* argv */)
{
  brick::utilities::TeeTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */

namespace {

  brick::utilities::TeeTest Currenttest;

}

#endif /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */
