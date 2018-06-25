/**
***************************************************************************
* @file brick/computerVision/test/getEuclideanDistanceTest.cc
*
* Source file defining tests for getEuclideanDistance().
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/getEuclideanDistance.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace computerVision {

    class GetEuclideanDistanceTest
      : public brick::test::TestFixture<GetEuclideanDistanceTest> {

    public:

      GetEuclideanDistanceTest();
      ~GetEuclideanDistanceTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testGetEuclideanDistance();

    private:

    }; // class GetEuclideanDistanceTest


    /* ============== Member Function Definititions ============== */

    GetEuclideanDistanceTest::
    GetEuclideanDistanceTest()
      : brick::test::TestFixture<GetEuclideanDistanceTest>("GetEuclideanDistanceTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testGetEuclideanDistance);
    }


    void
    GetEuclideanDistanceTest::
    testGetEuclideanDistance()
    {
      const size_t testImageRows = 240;
      const size_t testImageColumns = 320;
      std::vector<int> uCoords;
      std::vector<int> vCoords;

      uCoords.push_back(150);
      vCoords.push_back(100);

      uCoords.push_back(80);
      vCoords.push_back(130);

      uCoords.push_back(85);
      vCoords.push_back(170);

      uCoords.push_back(175);
      vCoords.push_back(165);

      // A simple test with only one "bright" pixel.
      Image<GRAY8> testImage(testImageRows, testImageColumns);
      testImage = common::UnsignedInt8(0);
      testImage(vCoords[0], uCoords[0]) = common::UnsignedInt8(1);

      numeric::Array2D<double> referenceDistances(testImage.rows(), testImage.columns());
      for(int row = 0; row < static_cast<int>(testImage.rows()); ++row) {
        for(int column = 0; column < static_cast<int>(testImage.columns());
            ++column) {
          double deltaU = column - uCoords[0];
          double deltaV = row - vCoords[0];
          double distance = std::sqrt(deltaU * deltaU + deltaV * deltaV);
          referenceDistances(row, column) = distance;
        }
      }

      numeric::Array2D<double> distances = getEuclideanDistance<double>(
        testImage, 10);

      BRICK_TEST_ASSERT(distances.rows() == referenceDistances.rows());
      BRICK_TEST_ASSERT(distances.columns() == referenceDistances.columns());
      for(size_t index0 = 0; index0 < distances.size(); ++index0) {
        BRICK_TEST_ASSERT(
          approximatelyEqual(distances[index0], referenceDistances[index0],
                             1.0E-6));
      }

      // A slightly harder test with additional bright pixels.
      for(size_t pointIndex = 1; pointIndex < uCoords.size(); ++pointIndex) {
        testImage(vCoords[pointIndex], uCoords[pointIndex]) = common::UnsignedInt8(1);
      }

      for(int row = 0; row < static_cast<int>(testImage.rows()); ++row) {
        for(int column = 0; column < static_cast<int>(testImage.columns());
            ++column) {
          double deltaU = column - uCoords[0];
          double deltaV = row - vCoords[0];
          double distance = std::sqrt(deltaU * deltaU + deltaV * deltaV);
          for(size_t pointIndex = 1; pointIndex < uCoords.size(); ++pointIndex) {
            deltaU = column - uCoords[pointIndex];
            deltaV = row - vCoords[pointIndex];
            distance =
              std::min(distance, std::sqrt(deltaU * deltaU + deltaV * deltaV));
          }
          referenceDistances(row, column) = distance;
        }
      }

      distances = getEuclideanDistance<double>(testImage, 10);

      BRICK_TEST_ASSERT(distances.rows() == referenceDistances.rows());
      BRICK_TEST_ASSERT(distances.columns() == referenceDistances.columns());
      for(size_t index0 = 0; index0 < distances.size(); ++index0) {
        BRICK_TEST_ASSERT(
          approximatelyEqual(distances[index0], referenceDistances[index0],
                             1.0E-6));
      }

    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::GetEuclideanDistanceTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::GetEuclideanDistanceTest currentTest;

}

#endif
