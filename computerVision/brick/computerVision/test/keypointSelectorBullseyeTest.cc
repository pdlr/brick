/**
***************************************************************************
* @file keypointSelectorBullseyeTest.cpp
*
* Source file defining tests for the KeypointSelectorBullseye class.
*
* Copyright (C) 2013 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <iomanip>
#include <iostream>

#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/keypointSelectorBullseye.hh>
#include <brick/computerVision/utilities.hh>
#include <brick/computerVision/test/testImages.hh>

#include <brick/test/functors.hh>
#include <brick/test/testFixture.hh>

#include <brick/utilities/timeUtilities.hh>

#include <brick/numeric/subArray2D.hh>

namespace brick {

  namespace computerVision {


    // This test is copied
    class KeypointSelectorBullseyeTest
      : public brick::test::TestFixture<KeypointSelectorBullseyeTest> {

    public:

      KeypointSelectorBullseyeTest();
      ~KeypointSelectorBullseyeTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testKeypointSelectorBullseye();
      void testExecutionTime();

    private:

      double m_defaultTolerance;

    }; // class KeypointSelectorBullseyeTest


    /* ============== Member Function Definititions ============== */

    KeypointSelectorBullseyeTest::
    KeypointSelectorBullseyeTest()
      : brick::test::TestFixture<KeypointSelectorBullseyeTest>("KeypointSelectorBullseyeTest"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testKeypointSelectorBullseye);
      // BRICK_TEST_REGISTER_MEMBER(testExecutionTime);
    }


    void
    KeypointSelectorBullseyeTest::
    testKeypointSelectorBullseye()
    {
      // We'll use several different input files.  The first is the
      // original test file, which contains a moderately tricky
      // bullseye.  The second is a copy of the first in which the
      // bullseye has been edited so that the centroid of the center
      // spot is at a precisely known location.  The third and fourth
      // are not-so-tricky images, but exercise the detector at
      // a variety of different scales.
      std::vector<std::string> inputFileNameVector;
      inputFileNameVector.push_back(getBullseyeFileNamePGM0());
      inputFileNameVector.push_back(getBullseyeFileNamePGM1());
      // inputFileNameVector.push_back(getBullseyeFileNamePGM2());
      // inputFileNameVector.push_back(getBullseyeFileNamePGM3());

      // How many bullseyes in each image.
      std::vector<std::size_t> bullseyeCounts;
      bullseyeCounts.push_back(1);
      bullseyeCounts.push_back(1);
      // bullseyeCounts.push_back(6);
      // bullseyeCounts.push_back(6);

      // How big the bullseyes are.
      std::vector<std::size_t> minRadiuses;
      minRadiuses.push_back(5);
      minRadiuses.push_back(5);
      // minRadiuses.push_back(3);
      // minRadiuses.push_back(30);

      std::vector<std::size_t> maxRadiuses;
      maxRadiuses.push_back(15);
      maxRadiuses.push_back(15);
      // maxRadiuses.push_back(15);
      // maxRadiuses.push_back(50);

      // These tolerances control how precisely the general position
      // output of the keypoint selector must match the nominal
      // position of the bullseye.  The tolerance is much tighter on
      // the second image because the centroid of the center bull of
      // the target is precisely known.  The tolerance is looser
      // on the 4th image because the bullseyes are large, making it
      // difficult to estimate their centroids by hand.
      std::vector<double> toleranceVector;
      toleranceVector.push_back(0.5);
      toleranceVector.push_back(1.0e-8);
      // toleranceVector.push_back(0.5);
      // toleranceVector.push_back(2.0);

      // Here is hand-calculated ground truth for each image.
      std::vector< std::vector<numeric::Vector2D<double> > >
        bullseyePositionVector;
      bullseyePositionVector.push_back(
        {numeric::Vector2D<double>(54.5, 59.5)}
      );
      bullseyePositionVector.push_back(
        {numeric::Vector2D<double>(54.5, 59.5)}
      );
      // bullseyePositionVector.push_back(
      //   {
      //     numeric::Vector2D<double>(108.5, 98.5),
      //     numeric::Vector2D<double>(168.7, 99.0),
      //     numeric::Vector2D<double>(228.5, 99.5),
      //     numeric::Vector2D<double>(107.5, 179.0),
      //     numeric::Vector2D<double>(167.6, 179.3),
      //     numeric::Vector2D<double>(227.7, 179.5)
      //   }
      // );
      // bullseyePositionVector.push_back(
      //   {
      //     numeric::Vector2D<double>(167.0, 165.0),
      //     numeric::Vector2D<double>(630.5, 173.0),
      //     numeric::Vector2D<double>(1085.0, 181.0),
      //     numeric::Vector2D<double>(161.0, 783.0),
      //     numeric::Vector2D<double>(624.0, 785.0),
      //     numeric::Vector2D<double>(1078.5, 787.0)
      //   }
      // );

      for(brick::common::UInt32 jj = 0; jj < inputFileNameVector.size();
	  ++jj) {

	// Load an image containing bullseye(s).
	std::string inputFileName = inputFileNameVector[jj];
	Image<GRAY8> inputImage = readPGM8(inputFileName);

	// Make sure the detector finds the targets.
	KeypointSelectorBullseye<double> selector(bullseyeCounts[jj],
                                                  maxRadiuses[jj],
                                                  minRadiuses[jj]);
	selector.setImage(inputImage);
	std::vector< KeypointBullseye<int> > keypoints =
	  selector.getKeypoints();


	BRICK_TEST_ASSERT(keypoints.size()
                          == bullseyePositionVector[jj].size());

        // Organize keypoints into raster order.
        std::sort(keypoints.begin(), keypoints.end(),
                  [](KeypointBullseye<int> const& xx,
                     KeypointBullseye<int> const& yy)
                  {return xx.column < yy.column;});
        std::sort(keypoints.begin(), keypoints.end(),
                  [](KeypointBullseye<int> const& xx,
                     KeypointBullseye<int> const& yy)
                  {return xx.row < yy.row;});

	for(brick::common::UInt32 ii = 0; ii < keypoints.size(); ++ii) {
          BRICK_TEST_ASSERT(
            brick::test::approximatelyEqual(
              keypoints[ii].row,
              static_cast<int>(bullseyePositionVector[jj][ii].y()),
              1));
          BRICK_TEST_ASSERT(
            brick::test::approximatelyEqual(
              keypoints[ii].column,
              static_cast<int>(bullseyePositionVector[jj][ii].x()),
              1));
        }

	// Make sure sub-pixel version is plausible.
	std::vector< KeypointBullseye<double> > keypointsGP =
	  selector.getKeypointsGeneralPosition();
	BRICK_TEST_ASSERT(keypointsGP.size() == keypoints.size());

        // Organize keypoints into raster order.
        std::sort(keypointsGP.begin(), keypointsGP.end(),
                  [](KeypointBullseye<double> const& xx,
                     KeypointBullseye<double> const& yy)
                  {return xx.column < yy.column;});
        std::sort(keypointsGP.begin(), keypointsGP.end(),
                  [](KeypointBullseye<double> const& xx,
                     KeypointBullseye<double> const& yy)
                  {return xx.row < yy.row;});
	for(brick::common::UInt32 ii = 0; ii < keypointsGP.size(); ++ii) {
          BRICK_TEST_ASSERT(
            brick::test::approximatelyEqual(
              keypointsGP[ii].row, bullseyePositionVector[jj][ii].y(),
              toleranceVector[jj]));
          BRICK_TEST_ASSERT(
            brick::test::approximatelyEqual(
              keypointsGP[ii].column, bullseyePositionVector[jj][ii].x(),
              toleranceVector[jj]));
	}
      }
    }

    void
    KeypointSelectorBullseyeTest::
    testExecutionTime()
    {
      // Load an image with a moderately tricky bullseye in it.
      std::string inputFileName = getBullseyeFileNamePGM0();
      Image<GRAY8> inputImage = readPGM8(inputFileName);

      brick::common::UInt32 const scale = 8;
      Image<GRAY8> bigImage(inputImage.rows() * scale,
                            inputImage.columns() * scale);
      for(brick::common::UInt32 ii = 0; ii < scale; ++ii) {
        numeric::subArray(bigImage,
                          numeric::Slice(ii * inputImage.rows(),
                                         (ii + 1) * inputImage.rows()),
                          numeric::Slice(ii * inputImage.columns(),
                                         (ii + 1) * inputImage.columns())) =
          numeric::subArray(inputImage);
      }

      writePGM8("big.pgm", bigImage);

      KeypointSelectorBullseye<double> selector(1, 15, 5);

      std::cout << "Start!" << std::endl;

      double t0 = brick::utilities::getCurrentTime();
      selector.setImage(bigImage);
      std::vector< KeypointBullseye<int> > keypoints = selector.getKeypoints();
      double t1 = brick::utilities::getCurrentTime();

      std::cout << "Stop!" << std::endl;
      std::cout << "Processed an image of size (" << bigImage.columns()
                << ", " << bigImage.rows() << ") in "
                << std::fixed << std::setprecision(5)
                << t1 - t0 << " seconds" << std::endl;
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main()
{
  brick::computerVision::KeypointSelectorBullseyeTest currentTest;
  bool result = currentTest.run();

  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::KeypointSelectorBullseyeTest currentTest;

}

#endif
