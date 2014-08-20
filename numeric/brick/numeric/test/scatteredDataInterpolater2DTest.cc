/**
***************************************************************************
* @file scatteredDataInterpolater2DTest.cpp
* 
* Source file defining ScatteredDataInterpolater2DTest class.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <stdint.h>

#include <limits>

#include <brick/common/functional.hh>
#include <brick/numeric/scatteredDataInterpolater2D.hh>
#include <brick/numeric/subArray2D.hh>
#include <brick/test/testFixture.hh>

using namespace brick::numeric;

namespace brick {

  namespace numeric {

    class ScatteredDataInterpolater2DTest : public brick::test::TestFixture<ScatteredDataInterpolater2DTest> {

    public:

      ScatteredDataInterpolater2DTest();
      ~ScatteredDataInterpolater2DTest() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testApproximate();

    private:

      double m_defaultTolerance;
    
    }; // class ScatteredDataInterpolater2DTest


    /* ============== Member Function Definititions ============== */

    ScatteredDataInterpolater2DTest::
    ScatteredDataInterpolater2DTest()
      : TestFixture<ScatteredDataInterpolater2DTest>("ScatteredDataInterpolater2DTest"),
        m_defaultTolerance(1.0E-7)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testApproximate);
    }


    void
    ScatteredDataInterpolater2DTest::
    testApproximate()
    {
      // Arbitrary, made up, scattered data to interpolate.
      Array2D<double> testData(
        "[[0.1, 0.7, 4.0],"
        " [5.5, 0.2, 6.0],"
        " [3.0, 3.0, 1.0],"
        " [1.2, 3.7, 0.4],"
        " [1.2, 3.6, 1.3],"
        " [4.2, 4.6, 6.0],"
        " [5.7, 2.5, 5.8],"
        " [3.1, 0.1, 2.2]]");

      Vector2D<double> corner0(0.0, 0.0);
      Vector2D<double> corner1(6.0, 5.0);

      // Break out the test data for easier access later.
      Array2D<double> sCoordArray = subArray(testData, Slice(), Slice(0, 1));
      Array2D<double> tCoordArray = subArray(testData, Slice(), Slice(1, 2));
      Array2D<double> valueArray = subArray(testData, Slice(), Slice(2, 3));
        
      // Figure out how many spline levels we need to perfectly match
      // each datapoint.
      double minimumPointSeparation = 0.1; // From inspection of testData.
      double span = std::max(corner1.getX() - corner0.getX(),
                             corner1.getY() - corner0.getY());
      uint32_t numberOfSamples = static_cast<uint32_t>(
        (span / minimumPointSeparation) + 0.5);

      // Interpolation uses B-Splines of order 3, so each sample has a
      // support of 4 control points.  Also there are additional
      // control points beyond the data range.
      uint32_t numberOfControlPoints = 4 * numberOfSamples + 3;

      // Resolution increases exponentially with number of levels.
      // Figure out which level results in the required number of
      // control points.
      uint32_t maxNumberOfLevels = static_cast<uint32_t>(
        std::log(numberOfControlPoints) / std::log(2) + 2);
      
      // Approximate the made up data using
      // ScatteredDataInterpolater2D instances of progressively
      // increasing resolution, verifying that the fidelity of the
      // approximation increases with each step.
      Array1D<double> currentResiduals(valueArray.size());
      Array1D<double> previousResiduals = valueArray.ravel().copy();
      for(uint32_t numberOfLevels = 1; numberOfLevels <= maxNumberOfLevels;
          ++numberOfLevels) {

        // Do the approximation.
        ScatteredDataInterpolater2D<double> scatteredDataInterpolater(
          numberOfLevels);
        scatteredDataInterpolater.approximate(
          sCoordArray.begin(), sCoordArray.end(), tCoordArray.begin(),
          valueArray.begin(), corner0, corner1);

        // Compute residuals so we can check that the approximation
        // improved (or at least didn't get worse) with this increase
        // in numberOfLevels.
        for(size_t index0 = 0; index0 < testData.rows(); ++index0) {
          double computedResult =
            scatteredDataInterpolater(testData(index0, 0), testData(index0, 1));
          currentResiduals[index0] = 
            common::absoluteValue(computedResult - testData(index0, 2));
        }

        double previousRMS = rms<double>(previousResiduals);
        double currentRMS = rms<double>(currentResiduals);
        BRICK_TEST_ASSERT(currentRMS < previousRMS + m_defaultTolerance);

        previousResiduals = currentResiduals.copy();
      }

      // Make sure we perfectly interpolated input data on the last
      // (highest resolution) try.
      double largestResidual = maximum<double>(currentResiduals);
      BRICK_TEST_ASSERT(
        approximatelyEqual(largestResidual, 0.0, this->m_defaultTolerance));

    } // testApproximate()
    
  } // namespace numeric

} // namespace brick


#if 0

int main(int /* argc */, char** /* argv */)
{
  brick::numeric::ScatteredDataInterpolater2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::ScatteredDataInterpolater2DTest currentTest;

}

#endif

