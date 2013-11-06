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

#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/kernels.hh>
#include <brick/computerVision/keypointSelectorBullseye.hh>
#include <brick/computerVision/utilities.hh>
#include <brick/computerVision/test/testImages.hh>

#include <brick/random/pseudoRandom.hh>
#include <brick/test/functors.hh>
#include <brick/test/testFixture.hh>

#include <brick/utilities/path.hh>
#include <brick/utilities/stringManipulation.hh>
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
      // Load an image with a moderately tricky bullseye in it.
      std::string inputFileName = getBullseyeFileNamePGM0();
      Image<GRAY8> inputImage = readPGM8(inputFileName);

      // Where we expect the keypoint detector to fire.
      // Note(xxx): must be a better way than hardcoding this.
      numeric::Index2D bullseyePosition(59, 54);
      
      // Make sure the detector finds the target.
      KeypointSelectorBullseye<double> selector(1, 15, 5);
      selector.setImage(inputImage);
      std::vector< KeypointBullseye<int> > keypoints = selector.getKeypoints();

      Image<GRAY8> flagImage(inputImage.rows(), inputImage.columns());
      flagImage = 0;
      for(brick::common::UInt32 ii = 0; ii < keypoints.size(); ++ii) {
        flagImage(keypoints[ii].row, keypoints[ii].column) =
          (ii + 1) * (255 / keypoints.size());
      }
      // writePGM8("flag.pgm", flagImage);

      BRICK_TEST_ASSERT(keypoints.size() == 1);
      BRICK_TEST_ASSERT(
        brick::test::approximatelyEqual(
          keypoints[0].row, bullseyePosition.getRow(), 1));
      BRICK_TEST_ASSERT(
        brick::test::approximatelyEqual(
          keypoints[0].column, bullseyePosition.getColumn(), 1));


      // Make sure sub-pixel version is plausible.
      std::vector< KeypointBullseye<double> > keypointsGP =
        selector.getKeypointsGeneralPosition();
      BRICK_TEST_ASSERT(keypointsGP.size() == 1);
      for(brick::common::UInt32 ii = 0; ii < keypointsGP.size(); ++ii) {
        BRICK_TEST_ASSERT(
          brick::test::approximatelyEqual(
            keypointsGP[ii].row,
            static_cast<double>(keypoints[ii].row) + 0.5,
            1.0));
        BRICK_TEST_ASSERT(
          brick::test::approximatelyEqual(
            keypointsGP[ii].column,
            static_cast<double>(keypoints[ii].column) + 0.5,
            1.0));
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
      
      // Make sure the detector finds the target.
      KeypointSelectorBullseye<double> selector(1, 15, 5);
      std::cout << "Start!" << std::endl;
      selector.setImage(bigImage);
      std::vector< KeypointBullseye<int> > keypoints = selector.getKeypoints();
      std::cout << "Stop!" << std::endl;
    }
    
  } // namespace computerVision

} // namespace brick


#if 0

int main()
{
  brick::computerVision::KeypointSelectorBullseyeTest currentTest;
  // currentTest.testKeypointSelectorBullseye();
  bool result = currentTest.run();

  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::KeypointSelectorBullseyeTest currentTest;

}

#endif
