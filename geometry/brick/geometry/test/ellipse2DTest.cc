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

// xxx
#include <iomanip>
#include <brick/utilities/imageIO.hh>

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
        m_defaultTolerance(1.0E-12)
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

        // xxx
        brick::numeric::Array1D<double> dRow(6);
        dRow[0] = xx * xx;
        dRow[1] = xx * root0;
        dRow[2] = root0 * root0;
        dRow[3] = xx;
        dRow[4] = root0;
        dRow[5] = 1.0;
        std::cout << std::fixed << std::setprecision(12)
                  << "dRow: " << dRow << "\n"
                  << "AD0: "
                  << brick::numeric::dot<double>(dRow, algebraicParameters)
                  << std::endl;
        
        samplePoints.push_back(
          brick::numeric::Vector2D<double>(xx, root0));
        samplePoints.push_back(
          brick::numeric::Vector2D<double>(xx, root1));
      }

      brick::numeric::Array1D<double> recoveredX;
      ellipse2D.estimate(samplePoints.begin(), samplePoints.end(),
                         recoveredX);

      // xxx
      brick::numeric::Array2D<unsigned char> graph(500, 500);
      graph = static_cast<unsigned char>(0);
      for(unsigned int ii = 0; ii < samplePoints.size(); ++ii) {
        unsigned int row = (samplePoints[ii].y() * 10 + 250.5);
        unsigned int column = (samplePoints[ii].x() * 10 + 250.5);
        if(row < graph.rows() && column < graph.columns()) {
          graph(row, column) = static_cast<unsigned char>(255);
        }
      }

#if 1
      for(double angle = 0.0; angle < 6.28; angle += (3.14 / 180.0)) {
        double ct = std::cos(angle);
        double st = std::sin(angle);
        brick::numeric::Vector2D<double> point = (ellipse2D.getOrigin()
                 + ct * ellipse2D.getSemimajorAxis()
                 + st * ellipse2D.getSemiminorAxis());
        unsigned int row = (point.y() * 10 + 250.5);
        unsigned int column = (point.x() * 10 + 250.5);
        if(row < graph.rows() && column < graph.columns()) {
          graph(row, column) = static_cast<unsigned char>(128);
        }
      }
#else
      std::cout << "recovered: " << recoveredX << std::endl;
      for(double xx = minimumX; xx < maximumX; xx += stepSize / 4) {
        double c0 = recoveredX[2];
        double c1 = recoveredX[1] * xx + recoveredX[4];
        double c2 = (recoveredX[0] * xx * xx
                     + recoveredX[3] * xx
                     + recoveredX[5]);
        double root0 = 0.0;
        double root1 = 0.0;
        bool isReal = brick::numeric::solveQuadratic(c0, c1, c2, root0, root1);
        if(!isReal) {
          // BRICK_THROW(brick::common::LogicException,
          //             "Ellipse2DTest::testEstimate()",
          //             "Valid xx range not valid for recoveredX.");
          std::cout << "(" << xx << ")" << std::flush;
        } else {

          brick::numeric::Vector2D<double> point0(xx, root0);
          unsigned int row = (point0.y() * 10 + 250.5);
          unsigned int column = (point0.x() * 10 + 250.5);
          if(row < graph.rows() && column < graph.columns()) {
            graph(row, column) = static_cast<unsigned char>(128);
          }
        
          brick::numeric::Vector2D<double> point1(xx, root1);
          row = (point1.y() * 10 + 250.5);
          column = (point1.x() * 10 + 250.5);
          if(row < graph.rows() && column < graph.columns()) {
            graph(row, column) = static_cast<unsigned char>(128);
          }
        }
      }
      std::cout << std::endl;
#endif      

      brick::utilities::writePGM("graph.pgm", graph.data(),
                                 graph.rows(), graph.columns());
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
