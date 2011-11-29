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
#include <brick/computerVision/test/testImages.hh>
#include <brick/test/testFixture.hh>

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
      Image<GRAY8> inputImage = readPGM8(getTestImageFileNamePGM0());
      KeypointSelectorFast selector;
      double time0 = utilities::getCurrentTime();
      selector.setImage(inputImage);
      double time1 = utilities::getCurrentTime();
      std::cout << "ET: " << time1 - time0 << " to search a "
                << inputImage.rows() << "x" << inputImage.columns()
                << " gray image." << std::endl;
      std::cout << "Threshold: " << selector.getThreshold() << std::endl;

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

int main(int argc, char** argv)
{
  brick::computerVision::KeypointSelectorFastTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::KeypointSelectorFastTest currentTest;

}

#endif
