/**
***************************************************************************
* @file bSpline2DTest.cpp
* 
* Source file defining BSpline2DTest class.
*
* Copyright (C) 2006-2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <limits>

#include <brick/common/functional.hh>
#include <brick/numeric/bSpline2D.hh>
#include <brick/numeric/subArray2D.hh>
#include <brick/test/testFixture.hh>

using namespace brick::numeric;

namespace brick {

  namespace numeric {

    class BSpline2DTest : public brick::test::TestFixture<BSpline2DTest> {

    public:

      BSpline2DTest();
      ~BSpline2DTest() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testApproximateScatteredData();

    private:

      double m_defaultTolerance;
    
    }; // class BSpline2DTest


    /* ============== Member Function Definititions ============== */

    BSpline2DTest::
    BSpline2DTest()
      : TestFixture<BSpline2DTest>("BSpline2DTest"),
        m_defaultTolerance(1.0E-7)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testApproximateScatteredData);
    }


    void
    BSpline2DTest::
    testApproximateScatteredData()
    {
      // Arbitrary, made up, scattered data to interpolate.
      Array2D<double> testData(
        "[[0.0, 0.7, 4.0],"
        " [5.5, 0.0, 6.0],"
        " [3.0, 3.0, 1.0],"
        " [1.2, 3.7, 0.4],"
        " [3.0, 4.5, 1.3],"
        " [4.2, 4.6, 6.0],"
        " [5.7, 2.5, 5.8],"
        " [3.1, 0.0, 2.2]]");


      // Approximate the made up data using a BSpline2D instance,
      // using a variety of control point densities.
      Array1D<double> previousResiduals(testData.rows());
      previousResiduals = std::numeric_limits<double>::max();
      Array1D<double> currentResiduals(testData.rows());
      for(size_t nn = 5; nn < 30; nn += 3) {

        BSpline2D<double> bSpline;
        bSpline.setNumberOfNodes(nn, nn + 2);
        Array2D<double> sCoordArray = subArray(testData, Slice(), Slice(0, 1));
        Array2D<double> tCoordArray = subArray(testData, Slice(), Slice(1, 2));
        Array2D<double> zCoordArray = subArray(testData, Slice(), Slice(2, 3));
        bSpline.approximateScatteredData(sCoordArray.begin(), sCoordArray.end(),
                                         tCoordArray.begin(),
                                         zCoordArray.begin());

        // Make sure approximation is better with this number of
        // control points than it was for the previous, coarser grid.
        // Note that this is not guaranteed to be true, but happens to
        // be true for the selected test data and control point
        // spacings.
        for(size_t index0 = 0; index0 < testData.rows(); ++index0) {
          double computedResult =
            bSpline(testData(index0, 0), testData(index0, 1));
          currentResiduals[index0] = 
            common::absoluteValue(computedResult - testData(index0, 2));

          std::cout << currentResiduals[index0] << " vs. " 
                    << previousResiduals[index0] << std::endl;
          BRICK_TEST_ASSERT(currentResiduals[index0]
                            <= previousResiduals[index0] + m_defaultTolerance);
        }

        // Get ready for the next control point spacing.
        previousResiduals = currentResiduals.copy();
      }

      // At the finest control point density, there should be no
      // compromise: all points reconstructed exactly.
      for(size_t index0 = 0; index0 < testData.rows(); ++index0) {
        BRICK_TEST_ASSERT(approximatelyEqual(currentResiduals[index0], 
                                             0.0, m_defaultTolerance));
      }

    } // testApproximateScatteredData.

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::BSpline2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::BSpline2DTest currentTest;

}

#endif

