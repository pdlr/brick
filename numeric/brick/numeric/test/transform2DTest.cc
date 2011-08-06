/**
***************************************************************************
* @file transform2DTest.cpp
* 
* Source file defining Transform2DTest class.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/transform2D.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class Transform2DTest : public test::TestFixture<Transform2DTest> {

    public:

      Transform2DTest();
      ~Transform2DTest() {};

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Tests of member functions.
      void testInvert();

    private:

    }; // class Transform2DTest


    /* ============== Member Function Definititions ============== */

    Transform2DTest::
    Transform2DTest()
      : TestFixture<Transform2DTest>("Transform2DTest")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testInvert);
    }


    void
    Transform2DTest::
    testInvert()
    {
      Transform2D<common::Float64> xf0(1.0, 2.0, 3.0,
                                       0.0, 4.0, 5.0,
                                       1.0, 0.0, 6.0);
      Transform2D<common::Float64> xf0Inverse = xf0.invert();
      Transform2D<common::Float64> ident = xf0 * xf0Inverse;
      double testEpsilon = 1.0e-12;
      for(size_t rowIndex = 0; rowIndex < 3; ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < 3; ++columnIndex) {
          if(rowIndex == columnIndex) {
            BRICK_TEST_ASSERT(
              approximatelyEqual(
                ident(rowIndex, columnIndex), 1.0, testEpsilon));
          } else {
            BRICK_TEST_ASSERT(
              approximatelyEqual(
                ident(rowIndex, columnIndex), 0.0, testEpsilon));
          }
        }
      }
    }

  } // namespace numeric
 
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::Transform2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Transform2DTest currentTest;

}

#endif
