/**
***************************************************************************
* @file bSpline2DTest.cpp
* 
* Source file defining BSpline2DTest class.
*
* Copyright (C) 2006-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

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

      // Parameters and interpolated data recovered from reference
      // BSpline implementation.  For now, we force ourselves to trust
      // these numbers.
      Array2D<double> referenceData(
        "[[ 0.39241379,  0.39241379,  3.62563104], "
        " [ 0.39241379,  0.39241379,  3.62563104], "
        " [ 0.39241379,  0.39241379,  3.62563104], "
        " [ 0.39241379,  0.39241379,  3.62563104], "
        " [ 1.76586207,  1.76586207,  1.30474756], "
        " [ 1.76586207,  1.76586207,  1.30474756], "
        " [ 1.76586207,  1.76586207,  1.30474756], "
        " [ 1.76586207,  1.76586207,  1.30474756], "
        " [ 3.13931034,  3.13931034,  2.14289096], "
        " [ 3.13931034,  3.13931034,  2.14289096], "
        " [ 3.13931034,  3.13931034,  2.14289096], "
        " [ 3.13931034,  3.13931034,  2.14289096], "
        " [ 4.51275862,  4.51275862,  5.7244739 ], "
        " [ 4.51275862,  4.51275862,  5.7244739 ], "
        " [ 4.51275862,  4.51275862,  5.7244739 ], "
        " [ 4.51275862,  4.51275862,  5.7244739 ]]");
    
      // Approximate the made up data using a BSpline2D instance.
      BSpline2D<double> bSpline;
      bSpline.setNumberOfNodes(6, 7);
      Array2D<double> sCoordArray = subArray(testData, Slice(), Slice(0, 1));
      Array2D<double> tCoordArray = subArray(testData, Slice(), Slice(1, 2));
      Array2D<double> zCoordArray = subArray(testData, Slice(), Slice(2, 3));
      bSpline.approximateScatteredData(sCoordArray.begin(), sCoordArray.end(),
                                       tCoordArray.begin(),
                                       zCoordArray.begin());

      // Compare the resulting approximation with data from our
      // reference implementation.
      for(size_t index0 = 0; index0 < referenceData.rows(); ++index0) {
        double computedResult =
          bSpline(referenceData(index0, 0), referenceData(index0, 1));
        BRICK_TEST_ASSERT(
          approximatelyEqual(computedResult, referenceData(index0, 2),
                             m_defaultTolerance));

      }
    }

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

