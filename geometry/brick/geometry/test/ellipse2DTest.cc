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

      // Given an algebraic parameterization of an ellipse and an X
      // coordinate, find the corresponding Y coordinates (if the X
      // coordinate is within the range of the ellipse, else return
      // false).
      bool
      solveAlgebraicEllipse(
        double xx,
        brick::numeric::Array1D<double> const& algebraicParameters,
        double& yy0,
        double& yy1);
      
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
      // Pick some ellipse parameters to recover.
      brick::numeric::Array1D<double> const algebraicParameters(
        "[3.0, -5.0, 7.0, -52.0, 23.0, 10.0]");

      // Make a wild assumption about where on the X axis, this
      // ellipse lives.  We'd like a nicer way to do this.
      double const minimumX = -200.0;
      double const maximumX = 200.0;
      Ellipse2D<double> ellipse2D;

      // Assuming our test polynomial is defined somewhere on x =
      // [minimumX, maximumX], find it more precisely by converting to
      // a single-variable quadratic in y, and looking for x values at
      // which it has real roots.
      double minimumXObserved = maximumX;
      double maximumXObserved = minimumX;
      for(double xx = -100.0; xx < 100.0; xx += 1.0) {
        double yy0 = 0.0;
        double yy1 = 0.0;
        if(this->solveAlgebraicEllipse(xx, algebraicParameters, yy0, yy1)) {
          minimumXObserved = std::min(minimumXObserved, xx);
          maximumXObserved = std::max(maximumXObserved, xx);
        }
      }

      if(minimumXObserved >= maximumX
         || maximumXObserved <= minimumX
         || minimumXObserved > maximumXObserved) {
        BRICK_THROW(brick::common::LogicException,
                    "Ellipse2DTest::testEstimate()",
                    "Test ellipse is not in the expected X range.");
      }

      double stepSize = (maximumXObserved - minimumXObserved) / 20;
      std::vector< brick::numeric::Vector2D<double> > samplePoints;
      for(double xx = minimumXObserved; xx < maximumXObserved; xx += stepSize) {
        double yy0 = 0.0;
        double yy1 = 0.0;
        if(!this->solveAlgebraicEllipse(xx, algebraicParameters, yy0, yy1)) {
          BRICK_THROW(brick::common::LogicException,
                      "Ellipse2DTest::testEstimate()",
                      "Valid xx range is suddenly invalid.");
        }

        samplePoints.push_back(
          brick::numeric::Vector2D<double>(xx, yy0));
        samplePoints.push_back(
          brick::numeric::Vector2D<double>(xx, yy1));
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

    
    bool
    Ellipse2DTest::
    solveAlgebraicEllipse(
      double xx,
      brick::numeric::Array1D<double> const& algebraicParameters,
      double& yy0,
      double& yy1)
    {
      double c0 = algebraicParameters[2];
      double c1 = algebraicParameters[1] * xx + algebraicParameters[4];
      double c2 = (algebraicParameters[0] * xx * xx
                   + algebraicParameters[3] * xx
                   + algebraicParameters[5]);
      return brick::numeric::solveQuadratic(c0, c1, c2, yy0, yy1);
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
