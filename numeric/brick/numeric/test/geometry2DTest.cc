/**
***************************************************************************
* @file geometry2DTest.cc
* 
* Source file defining tests for the geometry2DTest class.
*
* Copyright (C) 2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/geometry2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class Geometry2DTest : public brick::test::TestFixture<Geometry2DTest> {

    public:

      Geometry2DTest();
      ~Geometry2DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testBilaterate();

    private:

      const double m_defaultTolerance;
    
    }; // class Geometry2DTest


    /* ============== Member Function Definititions ============== */

    Geometry2DTest::
    Geometry2DTest()
      : brick::test::TestFixture<Geometry2DTest>("Geometry2DTest"),
        m_defaultTolerance(1.0E-6)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testBilaterate);
    }


    void
    Geometry2DTest::
    testBilaterate()
    {
      for(double x0 = -2.0; x0 < 2.5; x0 += 1.0) {
        for(double y0 = -2.0; y0 < 2.5; y0 += 1.0) {
          for(double x1 = -2.0; x1 < 2.5; x1 += 1.0) {
            for(double y1 = 7.0; y1 < 11.5; y1 += 1.0) {
              for(double x2 = 7.0; x2 < 11.5; x2 += 1.0) {
                for(double y2 = -2.0; y2 < 2.5; y2 += 1.0) {
                  Vector2D<double> center0(x0, y0);
                  Vector2D<double> center1(x1, y1);
                  Vector2D<double> target(x2, y2);
                  double r0 = magnitude<double>(target - center0);
                  double r1 = magnitude<double>(target - center1);
                  Vector2D<double> recoveredPoint0;
                  Vector2D<double> recoveredPoint1;
                
                  bool result = bilaterate<double>(
                    center0, center1, r0, r1, recoveredPoint0, recoveredPoint1);
                  if(!result) {
                    BRICK_TEST_ASSERT(false);
                  }

                  double distance0 = magnitude<double>(
                    recoveredPoint0 - target);
                  double distance1 = magnitude<double>(
                    recoveredPoint1 - target);
                  bool result2 = ((distance0 < m_defaultTolerance)
                                  || (distance1 < m_defaultTolerance));
                  if(!result2) {
                    BRICK_TEST_ASSERT(false);
                  }
                }
              }
            }
          }
        }
      }
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::Geometry2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Geometry2DTest currentTest;

}

#endif
