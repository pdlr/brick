/**
***************************************************************************
* @file keypointSelectorFastTest.cpp
*
* Source file defining tests for the KeypointSelectorFast class.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/keypointSelectorFast.hh>
#include <brick/computerVision/utilities.hh>

#include <brick/random/pseudoRandom.hh>
#include <brick/test/testFixture.hh>

#include <brick/utilities/path.hh>
#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace computerVision {

    
    // This test is copied 
    class KeypointSelectorFastTest
      : public brick::test::TestFixture<KeypointSelectorFastTest> {

    public:

      KeypointSelectorFastTest();
      ~KeypointSelectorFastTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testKeypointSelectorFast();

      // Legacy functions.
      void exerciseKeypointSelectorFast(std::string const& fileName);
        
    private:

      double m_defaultTolerance;
      
    }; // class KeypointSelectorFastTest


    /* ============== Member Function Definititions ============== */

    KeypointSelectorFastTest::
    KeypointSelectorFastTest()
      : TestFixture<KeypointSelectorFastTest>("KeypointSelectorFastTest"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testKeypointSelectorFast);
    }


    void
    KeypointSelectorFastTest::
    testKeypointSelectorFast()
    {
      // Pick a few arbitrary spots for keypoints.
      const unsigned int numberOfPoints = 4;
      numeric::Index2D points[numberOfPoints];
      points[0] = numeric::Index2D(10, 12);
      points[1] = numeric::Index2D(12, 57);
      points[2] = numeric::Index2D(23, 12);
      points[3] = numeric::Index2D(77, 67);

      // Create an image that's lots of below-threshold noise.
      Image<GRAY8> inputImage(100, 120);
      inputImage = 128;
      random::PseudoRandom prandom(0);
      for(unsigned int ii = 0; ii < inputImage.size(); ++ii) {
        inputImage[ii] += prandom.uniformInt(-5, 5);
      }

      // Now manually add four features of varying types.
      int bressenhamRows[16] = {-3, -3, -2, -1,  0,  1,  2,  3,
                                 3,  3,  2,  1,  0, -1, -2, -3};
      int bressenhamColumns[16] = {0,  1,  2,  3,  3,  3,  2,  1,
                                   0, -1, -2, -3, -3, -3, -2, -1};
      inputImage(points[0].getRow(), points[0].getColumn()) -= 40;
      inputImage(points[1].getRow(), points[1].getColumn()) += 40;
      inputImage(points[2].getRow(), points[2].getColumn()) -= 20;
      for(unsigned int ii = 14; ii < 14 + 16; ++ii) {
        unsigned int jj = ii % 16;
        inputImage(points[2].getRow() - bressenhamRows[jj],
                   points[2].getColumn() - bressenhamColumns[jj]) += 20;
      }
      inputImage(points[3].getRow(), points[3].getColumn()) += 20;
      for(unsigned int ii = 5; ii < 5 + 16; ++ii) {
        unsigned int jj = ii % 16;
        inputImage(points[3].getRow() - bressenhamRows[jj],
                   points[3].getColumn() - bressenhamColumns[jj]) -= 20;
      }

      // Make sure the detector picks up these four features in raster
      // order.
      KeypointSelectorFast selector;
      selector.setThreshold(29);
      selector.setImage(inputImage);
      std::vector<KeypointFast> keypoints = selector.getKeypoints();
      BRICK_TEST_ASSERT(keypoints.size() == numberOfPoints);
      for(unsigned int ii = 0; ii < numberOfPoints; ++ii) {
        BRICK_TEST_ASSERT(keypoints[ii].row == points[ii].getRow());
        BRICK_TEST_ASSERT(keypoints[ii].column == points[ii].getColumn());
        BRICK_TEST_ASSERT(keypoints[ii].isPositive == ((ii % 2) != 0));
        for(unsigned int jj = 0; jj < 16; ++jj) {
          common::UnsignedInt8 pixelValue = 
            inputImage(points[ii].getRow() + bressenhamRows[jj],
                       points[ii].getColumn() + bressenhamColumns[jj]);
          BRICK_TEST_ASSERT(keypoints[ii].featureVector[jj] == pixelValue);
        }
      }
    }


    void
    KeypointSelectorFastTest::
    exerciseKeypointSelectorFast(std::string const& fileName)
    {
      Image<GRAY8> inputImage;
      if(utilities::splitExtension(fileName).second == ".pgm") {
        inputImage = readPGM8(fileName);    
      } else {
        Image<RGB8> colorImage = readPPM8(fileName);
        inputImage = convertColorspace<GRAY8>(colorImage);
      }

      KeypointSelectorFast selector;
      double time0 = utilities::getCurrentTime();
      selector.setImage(inputImage);
      double time1 = utilities::getCurrentTime();
      std::cout << "ET: " << time1 - time0 << " to search a "
                << inputImage.rows() << "x" << inputImage.columns()
                << " gray image." << std::endl;
      std::cout << "Threshold: " << selector.getThreshold() << std::endl;
      time0 = utilities::getCurrentTime();
      selector.setImage(inputImage);
      time1 = utilities::getCurrentTime();
      std::cout << "2nd ET: " << time1 - time0 << std::endl;

      std::vector<KeypointFast> keypoints = selector.getKeypoints();
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

int main(/* int argc, char** argv */)
{
  brick::computerVision::KeypointSelectorFastTest currentTest;
  bool result = currentTest.run();

  currentTest.exerciseKeypointSelectorFast("testImagePGM0.pgm");
  
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::KeypointSelectorFastTest currentTest;
  
}

#endif
