/**
***************************************************************************
* @file bullseye2DTest.cpp
*
* Source file defining tests for the Bullseye2D class.
*
* Copyright (C) 2013 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/geometry/bullseye2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/solveQuadratic.hh>
#include <brick/test/testFixture.hh>

// xxx
#include <brick/geometry/ellipse2D.hh>
#include <brick/utilities/imageIO.hh>

namespace brick {

  namespace geometry {
    
    class Bullseye2DTest : public brick::test::TestFixture<Bullseye2DTest> {

    public:

      Bullseye2DTest();
      ~Bullseye2DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testEstimate();

      
    private:

      // xxx
      void drawGraphs();
      
      // Given an algebraic parameterization of an bullseye and an X
      // coordinate, find the corresponding Y coordinates (if the X
      // coordinate is within the range of the bullseye, else return
      // false).
      bool
      solveAlgebraicEllipse(
        double xx,
        brick::numeric::Array1D<double> const& algebraicParameters,
        double& yy0,
        double& yy1);
      
      const double m_defaultTolerance;
      
    }; // class Bullseye2DTest


    /* ============== Member Function Definititions ============== */

    Bullseye2DTest::
    Bullseye2DTest()
      : brick::test::TestFixture<Bullseye2DTest>("Bullseye2DTest"),
        m_defaultTolerance(1.0E-9)
    {
      BRICK_TEST_REGISTER_MEMBER(testEstimate);
    }


    void
    Bullseye2DTest::
    testEstimate()
    {
      this->drawGraphs();
      return;
      
      // Pick some bullseye parameters to recover.
      brick::numeric::Array1D<double> const algebraicParameters(
        "[3.0, -5.0, 7.0, -52.0, 23.0, 10.0]");

      // Make a wild assumption about where on the X axis, this
      // bullseye lives.  We'd like a nicer way to do this.
      double const minimumX = -200.0;
      double const maximumX = 200.0;
      Bullseye2D<double> bullseye2D;

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
                    "Bullseye2DTest::testEstimate()",
                    "Test bullseye is not in the expected X range.");
      }

      double stepSize = (maximumXObserved - minimumXObserved) / 20;
      std::vector< brick::numeric::Vector2D<double> > samplePoints;
      for(double xx = minimumXObserved; xx < maximumXObserved; xx += stepSize) {
        double yy0 = 0.0;
        double yy1 = 0.0;
        if(!this->solveAlgebraicEllipse(xx, algebraicParameters, yy0, yy1)) {
          BRICK_THROW(brick::common::LogicException,
                      "Bullseye2DTest::testEstimate()",
                      "Valid xx range is suddenly invalid.");
        }

        samplePoints.push_back(
          brick::numeric::Vector2D<double>(xx, yy0));
        samplePoints.push_back(
          brick::numeric::Vector2D<double>(xx, yy1));
      }

      // Recover the bullseye from the sample points.
      std::vector<unsigned int> counts;
      bullseye2D.estimate(samplePoints.begin(), samplePoints.end(),
                          counts.begin(), counts.end());

      // Verify that points on the bullseye match the parameterization
      // we started out with.
      for(double angle = 0.0; angle < 6.28; angle += (3.14 / 180.0)) {
        double ct = std::cos(angle);
        double st = std::sin(angle);
        brick::numeric::Vector2D<double> point = (bullseye2D.getOrigin()
                 + ct * bullseye2D.getSemimajorAxis(0)
                 + st * bullseye2D.getSemiminorAxis(0));

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

    
    void
    Bullseye2DTest::
    drawGraphs()
    {
      // Pick some bullseye parameters to recover.
      brick::numeric::Array2D<double> const algebraicParameters(
        "[[3.0, -5.0, 7.0, -52.0, 23.0, -100.0],"
        " [3.0, -5.0, 7.0, -52.0, 23.0, 10.0],"
        " [3.0, -5.0, 7.0, -52.0, 23.0, 200.0]]");
      brick::numeric::Array1D<brick::common::UInt8> colors(3);
      colors[0] = 255;
      colors[1] = 127;
      colors[2] = 63;

      brick::numeric::Array2D<unsigned char> graph(500, 500);
      graph = static_cast<unsigned char>(0);
      
      for(unsigned int rr = 0; rr < algebraicParameters.rows(); ++rr) {
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
          if(this->solveAlgebraicEllipse(
               xx, algebraicParameters.getRow(rr), yy0, yy1)) {
            minimumXObserved = std::min(minimumXObserved, xx);
            maximumXObserved = std::max(maximumXObserved, xx);
          }
        }

        if(minimumXObserved >= maximumX
           || maximumXObserved <= minimumX
           || minimumXObserved > maximumXObserved) {
          BRICK_THROW(brick::common::LogicException,
                      "Bullseye2DTest::drawGraphs()",
                      "Test ellipse is not in the expected X range.");
        }

        double stepSize = (maximumXObserved - minimumXObserved) / 20;
        std::vector< brick::numeric::Vector2D<double> > samplePoints;
        for(double xx = minimumXObserved; xx < maximumXObserved;
            xx += stepSize) {
          double yy0 = 0.0;
          double yy1 = 0.0;
          if(!this->solveAlgebraicEllipse(
               xx, algebraicParameters.getRow(rr), yy0, yy1)) {
            BRICK_THROW(brick::common::LogicException,
                        "Bullseye2DTest::drawGraphs()",
                        "Valid xx range is suddenly invalid.");
          }

          samplePoints.push_back(
            brick::numeric::Vector2D<double>(xx, yy0));
          samplePoints.push_back(
            brick::numeric::Vector2D<double>(xx, yy1));
        }

        for(unsigned int ii = 0; ii < samplePoints.size(); ++ii) {
          unsigned int row = (samplePoints[ii].y() * 10 + 250.5);
          unsigned int column = (samplePoints[ii].x() * 10 + 250.5);
          if(row < graph.rows() && column < graph.columns()) {
            graph(row, column) = colors[rr];
          }
        }
      }

      brick::utilities::writePGM("graph.pgm", graph.data(),
                                 graph.rows(), graph.columns());
    }


    bool
    Bullseye2DTest::
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
  brick::geometry::Bullseye2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Bullseye2DTest currentTest;

}

#endif
