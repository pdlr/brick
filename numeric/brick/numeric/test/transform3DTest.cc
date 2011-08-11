/**
***************************************************************************
* @file transform3DTest.cpp
* 
* Source file defining Transform3DTest class.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/transform3D.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {
    
    class Transform3DTest : public test::TestFixture<Transform3DTest> {

    public:

      Transform3DTest();
      ~Transform3DTest() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testInvert();

    private:

    }; // class Transform3DTest


    /* ============== Member Function Definititions ============== */

    Transform3DTest::
    Transform3DTest()
      : TestFixture<Transform3DTest>("Transform3DTest")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testInvert);
    }


    void
    Transform3DTest::
    testInvert()
    {
      Transform3D<common::Float64> xf0(1.0, 2.0, 3.0, 4.0,
                                       0.0, 3.2, -1.4, 11.0,
                                       -5.0, 0.0, 4.0, 6.0,
                                       2.0, -1.0, 2.0, 0.5);
      Transform3D<common::Float64> xf0Inverse = xf0.invert();
      Transform3D<common::Float64> ident = xf0 * xf0Inverse;
      double testEpsilon = 1.0e-12;
      for(size_t rowIndex = 0; rowIndex < 4; ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < 4; ++columnIndex) {
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
  brick::numericTransform3DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Transform3DTest currentTest;

}

#endif
