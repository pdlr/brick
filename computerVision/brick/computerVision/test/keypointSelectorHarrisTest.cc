/**
***************************************************************************
* @file keypointSelectorHarrisTest.cpp
*
* Source file defining tests for the KeypointSelectorHarris class.
*
* Copyright (C) 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/keypointSelectorHarris.hh>
#include <brick/computerVision/utilities.hh>

#include <brick/random/pseudoRandom.hh>
#include <brick/test/functors.hh>
#include <brick/test/testFixture.hh>

#include <brick/utilities/path.hh>
#include <brick/utilities/stringManipulation.hh>
#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace computerVision {

    
    // This test is copied 
    class KeypointSelectorHarrisTest
      : public brick::test::TestFixture<KeypointSelectorHarrisTest> {

    public:

      KeypointSelectorHarrisTest();
      ~KeypointSelectorHarrisTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testKeypointSelectorHarris();

      // Legacy functions.
      void exerciseKeypointSelectorHarris(std::string const& fileName,
                                          double threshold = 0.0);
        
    private:

      double m_defaultTolerance;
      
    }; // class KeypointSelectorHarrisTest


    /* ============== Member Function Definititions ============== */

    KeypointSelectorHarrisTest::
    KeypointSelectorHarrisTest()
      : TestFixture<KeypointSelectorHarrisTest>("KeypointSelectorHarrisTest"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testKeypointSelectorHarris);
    }


    void
    KeypointSelectorHarrisTest::
    testKeypointSelectorHarris()
    {
      // Create an image with well-defined corners in it.
      Image<GRAY8> inputImage(100, 120);
      inputImage = brick::common::UInt8(60);
      inputImage.getROI(numeric::Index2D(20, 30), numeric::Index2D(60, 70)) =
        common::UInt8(128);

      // Write down what our groundtruth ought to be for these four points.
      unsigned int const numberOfPoints = 4;
      std::vector<numeric::Index2D> points(numberOfPoints);
      points[0].setValue(60, 70);
      points[1].setValue(60, 30);
      points[2].setValue(20, 70);
      points[3].setValue(20, 30);

      // Where we expect the keypoint detector to fire.
      std::vector<numeric::Index2D> referencePoints(numberOfPoints);
      referencePoints[0].setValue(58, 68);
      referencePoints[1].setValue(58, 31);
      referencePoints[2].setValue(21, 68);
      referencePoints[3].setValue(21, 31);
      
      // Make sure the detector picks up these four corners.
      KeypointSelectorHarris<double> selector;
      selector.setImage(inputImage);
      std::vector< KeypointHarris<int> > keypoints = selector.getKeypoints();
      BRICK_TEST_ASSERT(keypoints.size() == numberOfPoints);
      for(unsigned int ii = 0; ii < numberOfPoints; ++ii) {
        BRICK_TEST_ASSERT(keypoints[ii].row
                          == referencePoints[ii].getRow());
        BRICK_TEST_ASSERT(keypoints[ii].column
                          == referencePoints[ii].getColumn());
      }


      // Make sure sub-pixel versions are plausible.
      std::vector< KeypointHarris<double> > keypointsGP =
        selector.getKeypointsGeneralPosition();
      BRICK_TEST_ASSERT(keypoints.size() == numberOfPoints);
      for(unsigned int ii = 0; ii < numberOfPoints; ++ii) {
        BRICK_TEST_ASSERT(
          brick::test::approximatelyEqual(
            keypointsGP[ii].row, (double)referencePoints[ii].getRow(),
            0.5));
        BRICK_TEST_ASSERT(
          brick::test::approximatelyEqual(
            keypointsGP[ii].column, (double)referencePoints[ii].getColumn(),
            0.5));
      }
    }


    void
    KeypointSelectorHarrisTest::
    exerciseKeypointSelectorHarris(std::string const& fileName,
                                   double threshold)
    {
      Image<GRAY8> inputImage;
      if(utilities::splitExtension(fileName).second == ".pgm") {
        inputImage = readPGM8(fileName);    
      } else {
        Image<RGB8> colorImage = readPPM8(fileName);
        inputImage = convertColorspace<GRAY8>(colorImage);
      }

      KeypointSelectorHarris<double> selector;
      std::vector< KeypointHarris<common::Int32> > keypoints;

      keypoints.clear();
      double time0 = utilities::getCurrentTime();
      selector.setImage(inputImage);
      // keypoints = selector.getKeypoints();
      selector.getKeypoints(std::back_inserter(keypoints), threshold);
      double time1 = utilities::getCurrentTime();
      std::cout << "ET: " << time1 - time0 << " to search a "
                << inputImage.rows() << "x" << inputImage.columns()
                << " gray image." << std::endl;

      keypoints.clear();
      time0 = utilities::getCurrentTime();
      selector.setImage(inputImage);
      // keypoints = selector.getKeypoints(threshold);
      selector.getKeypoints(std::back_inserter(keypoints), threshold);
      time1 = utilities::getCurrentTime();
      std::cout << "2nd ET: " << time1 - time0 << std::endl;

      std::cout << "Num keypoints: " << keypoints.size() << std::endl;
      for(unsigned int ii = 0; ii < keypoints.size(); ++ii) {
        brick::numeric::Index2D
          keyPoint(keypoints[ii].row, keypoints[ii].column);
        inputImage(keyPoint.getRow() - 2, keyPoint.getColumn() - 1) = 255;
        inputImage(keyPoint.getRow() - 2, keyPoint.getColumn()) = 255;
        inputImage(keyPoint.getRow() - 2, keyPoint.getColumn() + 1) = 255;
        inputImage(keyPoint.getRow() - 1, keyPoint.getColumn() - 2) = 255;
        inputImage(keyPoint.getRow() - 1, keyPoint.getColumn() + 2) = 255;
        inputImage(keyPoint.getRow(), keyPoint.getColumn() - 2) = 255;
        inputImage(keyPoint.getRow(), keyPoint.getColumn() + 2) = 255;
        inputImage(keyPoint.getRow() + 1, keyPoint.getColumn() - 2) = 255;
        inputImage(keyPoint.getRow() + 1, keyPoint.getColumn() + 2) = 255;
        inputImage(keyPoint.getRow() + 2, keyPoint.getColumn() - 1) = 255;
        inputImage(keyPoint.getRow() + 2, keyPoint.getColumn()) = 255;
        inputImage(keyPoint.getRow() + 2, keyPoint.getColumn() + 1) = 255;
      }

      writePGM("keypointImage.pgm", inputImage.data(),
               inputImage.rows(), inputImage.columns(), false);
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::KeypointSelectorHarrisTest currentTest;
  bool result = currentTest.run();

  double threshold = 0.0;
  if(argc > 1) {
    threshold = brick::utilities::convertString<double>(argv[1]);
  }
  // currentTest.exerciseKeypointSelectorHarris("testImagePGM0.pgm");
  currentTest.exerciseKeypointSelectorHarris("rimg.pgm", threshold);

  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::KeypointSelectorHarrisTest currentTest;

}

#endif
