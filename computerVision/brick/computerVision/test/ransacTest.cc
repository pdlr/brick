/**
***************************************************************************
* @file ransacTest.cpp
*
* Source file defining tests for RANSAC helper functions.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <algorithm>
#include <functional>
#include <brick/computerVision/ransac.hh>
#include <brick/computerVision/ransacClassInterface.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/random/pseudoRandom.hh>
#include <brick/test/testFixture.hh>

namespace num = brick::numeric;
namespace rndm = brick::random;

namespace brick {

  namespace computerVision {

    
    // This test is copied 
    class RansacTest
      : public brick::test::TestFixture<RansacTest> {

    public:

      RansacTest();
      ~RansacTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testRansac();

    private:

      double m_defaultTolerance;
      
    }; // class RansacTest


    /* ===== Declarations for a simple Ransac Problem ====== */

    // The first template argument says that we're going to estimate
    // the line based on a collection of Vector2D instances.  The
    // second template argument says that we're going to represent a
    // line using a pair of doubles (slope, intercept).
    class LineFittingProblem
      : public RansacProblem< num::Vector2D<double>, std::pair<double, double> >
    {
    public:

      // When constructing the LineFittingProblem instance, we have to
      // pass in a sequence of Vector2D instances so that the RANSAC
      // implementation can later ask our LineFittingProblem instance
      // for randomly selected sets of samples.  The first argument of
      // the parent class constructor, 2, indicates that two samples
      // (Vector2D instances) are needed to estimate a line.
      template <class IterType>
      LineFittingProblem(IterType beginIter, IterType endIter)
        : RansacProblem< num::Vector2D<double>, std::pair<double, double> >(
            2, beginIter, endIter) {}


      // SampleSequenceType is typedef'd inside RansacProblem.  It's
      // just a pair of iterators that define a sequence of Vector2D
      // instances.
      std::pair<double, double>
      estimateModel(SampleSequenceType const& sampleSequence) {
        // To estimate slope and intercept from a bunch of (x, y) pairs, we
        // have:
        //
        // @code
        //   |y_0|       |x_0|
        //   |y_1|       |x_1|
        //   |...| = m * |...| + b
        //   |y_n|       |x_n|
        // @endcode
        //
        // Rearranging gives:
        //
        // @code
        //   |dot(X,X), sum(X)|   |m|   |dot(X,Y)|
        //   |                | * | | = |        |
        //   | sum(X),    n   |   |b|   | sum(Y) |
        // @endocode
        //
        // where X = [x_0, x_1, ..., x_n]^T and Y = [y_0, y_1, ...,
        // y_n]^T, which we solve by cofactor inversion:
        //
        // @code
        //   |m}             1             |   n   , -sum(X) |   |dot(X,Y)|
        //   | | = --------------------- * |                 | * |        |
        //   |b|   n*dot(X,X) - sum(X)^2   |-sum(X), dot(X,X)|   | sum(Y) |
        // @endocode

        double dotXX = 0.0;
        double dotXY = 0.0;
        double sumX = 0.0;
        double sumY = 0.0;
        size_t nn = 0;

        // We have to copy sampleSequence so that we can increment the
        // "begin iterator" part of it.
        SampleSequenceType mutableSequence = sampleSequence;
        
        while(mutableSequence.first != mutableSequence.second) {
          num::Vector2D<double> const& xy_n = *mutableSequence.first;
          dotXX += xy_n.x() * xy_n.x();
          dotXY += xy_n.x() * xy_n.y();
          sumX += xy_n.x();
          sumY += xy_n.y();
          ++nn;

          ++mutableSequence.first;
        }
        double determinant = nn * dotXX - sumX * sumX;
        if(determinant == 0.0) {
          return std::make_pair(std::numeric_limits<double>::max(),
                                std::numeric_limits<double>::max());
        }
        double slope = (nn * dotXY - sumX * sumY) / determinant;
        double intercept = (dotXX * sumY - dotXY * sumX) / determinant;
        return std::make_pair(slope, intercept);
      }


      // Given a model (a line), compute the residual for each element
      // of the input sampleSequence.
      template <class IterType>
      void
      computeError(std::pair<double, double> const& model,
                   SampleSequenceType const& sampleSequence,
                   IterType outputIter) {
        SampleSequenceType mutableSequence = sampleSequence;
        while(mutableSequence.first != mutableSequence.second) {
          
          // We're looking for the point that minimizes squared error.
          // We parameterize this point by its x coord, and minimize
          // over x:
          //
          // @code
          // x_best = argmin_x(
          //   (sample.x() - x)**2 + (sample.y() - (intercept + x * slope))**2)
          // @endcode
          //
          // Since this is quadratic and concave up, we can simply find
          // the x that makes its first derivative equal to zero.
          //
          // @code
          // 0 = ((-2 * (sample.x() - x)
          //       - 2 * slope * (sample.y() - (intercept + x * slope)))
          // @endcode
          //
          // Solving for x, this gives
          //
          // @code
          // x = (slope * (sample.y() - intercept) + sample.x())
          //     / (slope**2 + 1)
          // @endcode

          num::Vector2D<double> sample = *mutableSequence.first;
          double slope = model.first;
          double intercept = model.second;
          double xBest = ((slope * (sample.y() - intercept) + sample.x())
                          / (slope * slope + 1));
          double yBest = intercept + slope * xBest;
          *outputIter = num::magnitude<double>(sample - num::Vector2D<double>(xBest, yBest));

          ++outputIter;
          ++mutableSequence.first;
        }
      }

      // How close must a point be to the line in order to not be
      // considered an outlier.
      double 
      getNaiveErrorThreshold() {return 0.5;}
    
    };


    /* ============== Member Function Definititions ============== */

    RansacTest::
    RansacTest()
      : brick::test::TestFixture<RansacTest>("RansacTest"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testRansac);
    }


    void
    RansacTest::
    testRansac()
    {
      std::vector< num::Vector2D<double> > sampleVector;
      sampleVector.push_back(num::Vector2D<double>(0.0, 0.0));
      sampleVector.push_back(num::Vector2D<double>(1.0, 1.0));
      sampleVector.push_back(num::Vector2D<double>(2.0, 2.0));
      sampleVector.push_back(num::Vector2D<double>(3.0, 2.0));
      sampleVector.push_back(num::Vector2D<double>(3.0, 3.0));
      sampleVector.push_back(num::Vector2D<double>(4.0, 4.0));
      sampleVector.push_back(num::Vector2D<double>(10.0, 2.0));

      // Require 5 inliers to terminate.
      size_t minConsensusSetSize = 5;

      // Only fail this test every 10 billion iterations or so.
      // Note(xxx): Make a way to say "iterate forever if necessary."
      double requiredConfidence = 1.0 - 1.0E-10;

      // Five out of seven sample points are inliers.
      double inlierProbability = 5.0 / 7.0;

      // Make the ransac instance to test.
      LineFittingProblem lineFittingProblem(
        sampleVector.begin(), sampleVector.end());
      Ransac<LineFittingProblem> ransac(
        lineFittingProblem, minConsensusSetSize, requiredConfidence,
        inlierProbability);

      std::pair<double, double> slope_intercept = ransac.getResult();
      BRICK_TEST_ASSERT(approximatelyEqual(slope_intercept.first, 1.0,
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(slope_intercept.second, 0.0,
                                         m_defaultTolerance));
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::RansacTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::RansacTest currentTest;

}

#endif
