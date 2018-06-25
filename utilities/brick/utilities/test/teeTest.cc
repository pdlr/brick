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
      Tee<char> tee;
      std::ostringstream stream0;
      std::ostringstream stream1;
      std::ostringstream stream2;

      std::size_t ii = 0;
      tee.add(stream0);
      tee << "Test Number " << ++ii << " " << std::endl;
      tee.add(stream1);
      tee << "Test Number " << ++ii << " " << std::endl;
      tee.add(stream2);
      tee << "Test Number " << ++ii << " " << std::endl;

      BRICK_TEST_ASSERT(
        stream0.str() == "Test Number 1 \nTest Number 2 \nTest Number 3 \n");
      BRICK_TEST_ASSERT(
        stream1.str() == "Test Number 2 \nTest Number 3 \n");
      BRICK_TEST_ASSERT(
        stream2.str() == "Test Number 3 \n");
    }

  } // namespace utilities

} // namespace brick


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
