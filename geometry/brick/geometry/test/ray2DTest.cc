/**
***************************************************************************
* @file ray2DTest.cpp
*
* Source file defining tests for Ray2D class.
*
* Copyright (C) 2009 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <cmath>
#include <brick/geometry/ray2D.hh>
#include <brick/test/testFixture.hh>

namespace num = brick::numeric;

namespace brick {

  namespace geometry {
    
    class Ray2DTest : public TestFixture<Ray2DTest> {

    public:

      Ray2DTest();
      ~Ray2DTest() {}

      void setUp(const std::string& testName) {}
      void tearDown(const std::string& testName) {}

      // Tests.
      void testConstructor__double__double__double();

    private:

      const double m_defaultTolerance;
      
    }; // class Ray2DTest


    /* ============== Member Function Definititions ============== */

    Ray2DTest::
    Ray2DTest()
      : TestFixture<Ray2DTest>("Ray2DTest"),
        m_defaultTolerance(1.0E-12)
    {
      BRICK_TEST_REGISTER_MEMBER(testConstructor__double__double__double);
    }


    void
    Ray2DTest::
    testConstructor__double__double__double()
    {
      double aa = 1.0;
      double bb = -2.0;
      double cc = 3.0;

      Ray2D ray2D(aa, bb, cc);
      for(size_t ii = 0; ii < 10; ++ii) {
        num::Vector2D testPoint =
          ray2D.getOrigin() + ii * ray2D.getDirectionVector();
        double residual = aa * testPoint.x() + bb * testPoint.y() + cc;
        BRICK_TEST_ASSERT(std::fabs(residual) < m_defaultTolerance);
      }
    }

  } // namespace geometry

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::geometry::Ray2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Ray2DTest currentTest;

}

#endif
