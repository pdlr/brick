/**
***************************************************************************
* @file keypointSelectorBullseyeTest.cpp
*
* Source file defining tests for the KeypointSelectorBullseye class.
*
* Copyright (C) 2013 David LaRose, dlr@cs.cmu.edu
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
      // We'll use two different input files.  The first is the
      // original test file, which contains a moderately tricky
      // bullseye.  The second is a copy of the first in which the
      // bullseye has been edited so that the centroid of the center
      // spot is at a precisely known location.
      std::vector<std::string> inputFileNameVector;
      inputFileNameVector.push_back(getBullseyeFileNamePGM0());
      inputFileNameVector.push_back(getBullseyeFileNamePGM1());

      // These tolerances control how precisely the general position
      // output of the keypoint selector must match the nominal
      // position of the bullseye.  The tolerance is much tighter on
      // the second image because the centroid of the center bull of
      // the target is precisely known.
      std::vector<double> toleranceVector;
      toleranceVector.push_back(0.5);
      toleranceVector.push_back(1.0e-8);

      // For now, both images have the same nominal bullseye position.
      std::vector< numeric::Vector2D<double> > bullseyePositionVector;
      bullseyePositionVector.push_back(numeric::Vector2D<double>(54.5, 59.5));
      bullseyePositionVector.push_back(numeric::Vector2D<double>(54.5, 59.5));
      
      for(brick::common::UInt32 jj = 0; jj < inputFileNameVector.size();
	  ++jj) {

	// Load an image with a moderately tricky bullseye in it.
	std::string inputFileName = inputFileNameVector[jj];
	Image<GRAY8> inputImage = readPGM8(inputFileName);

	// Make sure the detector finds the target.
	KeypointSelectorBullseye<double> selector(1, 15, 5);
	selector.setImage(inputImage);
	std::vector< KeypointBullseye<int> > keypoints =
	  selector.getKeypoints();


	BRICK_TEST_ASSERT(keypoints.size() == 1);
	BRICK_TEST_ASSERT(
          brick::test::approximatelyEqual(
            keypoints[0].row,
	    static_cast<int>(bullseyePositionVector[jj].y()),
	    1));
	BRICK_TEST_ASSERT(
          brick::test::approximatelyEqual(
            keypoints[0].column,
	    static_cast<int>(bullseyePositionVector[jj].x()),
	    1));

	// Make sure sub-pixel version is plausible.
	std::vector< KeypointBullseye<double> > keypointsGP =
	  selector.getKeypointsGeneralPosition();
	BRICK_TEST_ASSERT(keypointsGP.size() == 1);

	for(brick::common::UInt32 ii = 0; ii < keypointsGP.size(); ++ii) {
          BRICK_TEST_ASSERT(
            brick::test::approximatelyEqual(
              keypointsGP[0].row, bullseyePositionVector[jj].y(),
              toleranceVector[jj]));
          BRICK_TEST_ASSERT(
            brick::test::approximatelyEqual(
              keypointsGP[0].column, bullseyePositionVector[jj].x(),
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
