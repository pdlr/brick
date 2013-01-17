/**
***************************************************************************
* @file ellipse2DTest.cpp
*
* Source file defining tests for the Ellipse2D class.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/geometry/ellipse2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/solveQuadratic.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace geometry {
    
    class Ellipse2DTest : public brick::test::TestFixture<Ellipse2DTest> {

    public:

      Ellipse2DTest();
      ~Ellipse2DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testEstimate();

    private:

      const double m_defaultTolerance;
      
    }; // class Ellipse2DTest


    /* ============== Member Function Definititions ============== */

    Ellipse2DTest::
    Ellipse2DTest()
      : brick::test::TestFixture<Ellipse2DTest>("Ellipse2DTest"),
        m_defaultTolerance(1.0E-9)
    {
      BRICK_TEST_REGISTER_MEMBER(testEstimate);
    }


    void
    Ellipse2DTest::
    testEstimate()
    {
      brick::numeric::Array1D<double> algebraicParameters(
        "[3.0, -5.0, 7.0, -52.0, 23.0, 10.0]");
      Ellipse2D<double> ellipse2D;

      // Assuming our test polynomial is defined somewhere on x =
      // [-100, 100], find the valid range by converting to a
      // single-variable quadratic in y, and looking for x values at
      // which it has real roots.
      double minimumX = 200.0;
      double maximumX = -200.0;
      for(double xx = -100.0; xx < 100.0; xx += 1.0) {
        double c0 = algebraicParameters[2];
        double c1 = algebraicParameters[1] * xx + algebraicParameters[4];
        double c2 = (algebraicParameters[0] * xx * xx
                     + algebraicParameters[3] * xx
                     + algebraicParameters[5]);
        double root0 = 0.0;
        double root1 = 0.0;
        bool isReal = brick::numeric::solveQuadratic(c0, c1, c2, root0, root1);
        if(isReal) {
          minimumX = std::min(minimumX, xx);
          maximumX = std::max(maximumX, xx);
        }
      }

      if(minimumX > 100 || maximumX < -100 || minimumX > maximumX) {
        BRICK_THROW(brick::common::LogicException,
                    "Ellipse2DTest::testEstimate()",
                    "Test ellipse is not inside X range of [-100, 100].");
      }

      double stepSize = (maximumX - minimumX) / 20;
      std::vector< brick::numeric::Vector2D<double> > samplePoints;
      for(double xx = minimumX; xx < maximumX; xx += stepSize) {
        double c0 = algebraicParameters[2];
        double c1 = algebraicParameters[1] * xx + algebraicParameters[4];
        double c2 = (algebraicParameters[0] * xx * xx
                     + algebraicParameters[3] * xx
                     + algebraicParameters[5]);
        double root0 = 0.0;
        double root1 = 0.0;
        bool isReal = brick::numeric::solveQuadratic(c0, c1, c2, root0, root1);
        if(!isReal) {
          BRICK_THROW(brick::common::LogicException,
                      "Ellipse2DTest::testEstimate()",
                      "Valid xx range appears invalid.");
        }

        samplePoints.push_back(
          brick::numeric::Vector2D<double>(xx, root0));
        samplePoints.push_back(
          brick::numeric::Vector2D<double>(xx, root1));
      }

      // Recover the ellipse from the sample points.
      ellipse2D.estimate(samplePoints.begin(), samplePoints.end());

      // Verify that points on the ellipse match the parameterization
      // we started out with.
      for(double angle = 0.0; angle < 6.28; angle += (3.14 / 180.0)) {
        double ct = std::cos(angle);
        double st = std::sin(angle);
        brick::numeric::Vector2D<double> point = (ellipse2D.getOrigin()
                 + ct * ellipse2D.getSemimajorAxis()
                 + st * ellipse2D.getSemiminorAxis());

        double algebraicDistance = (
          algebraicParameters[0] * point.x() * point.x()
          + algebraicParameters[1] * point.x() * point.y()
          + algebraicParameters[2] * point.y() * point.y()
          + algebraicParameters[3] * point.x()
          + algebraicParameters[4] * point.y()
          + algebraicParameters[5]);

        BRICK_TEST_ASSERT(algebraicDistance < this->m_defaultTolerance);
      }
    }

  } // namespace geometry

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::geometry::Ellipse2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Ellipse2DTest currentTest;

}

#endif
