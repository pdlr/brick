/**
***************************************************************************
* @file keypointSelectorLepetitTest.cpp
*
* Source file defining tests for the KeypointSelectorLepetit class.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/keypointSelectorLepetit.hh>
#include <brick/computerVision/test/testImages.hh>
#include <brick/test/testFixture.hh>

#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace computerVision {

    
    // This test is copied 
    class KeypointSelectorLepetitTest
      : public brick::test::TestFixture<KeypointSelectorLepetitTest> {

    public:

      KeypointSelectorLepetitTest();
      ~KeypointSelectorLepetitTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testKeypointSelectorLepetit();

    private:

      double m_defaultTolerance;
      
    }; // class KeypointSelectorLepetitTest


    /* ============== Member Function Definititions ============== */

    KeypointSelectorLepetitTest::
    KeypointSelectorLepetitTest()
      : brick::test::TestFixture<KeypointSelectorLepetitTest>("KeypointSelectorLepetitTest"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testKeypointSelectorLepetit);
    }


    void
    KeypointSelectorLepetitTest::
    testKeypointSelectorLepetit()
    {
      Image<GRAY8> inputImage = readPGM8(getTestImageFileNamePGM0());
      double time0 = utilities::getCurrentTime();
      KeypointSelectorLepetit selector;
      selector.setImage(inputImage);
      double time1 = utilities::getCurrentTime();
      std::cout << "ET: " << time1 - time0 << " to search a "
                << inputImage.rows() << "x" << inputImage.columns()
                << " gray image." << std::endl;

      std::cout << "Pixel similarity threshold: "
                << selector.getPixelSimilarityThreshold() << "\n"
                << "Laplacian magnitude threshold: "
                << selector.getLaplacianMagnitudeThreshold() << std::endl;

      std::vector<brick::numeric::Index2D> keyPoints = selector.getKeypoints();
      for(unsigned int ii = 0; ii < keyPoints.size(); ++ii) {
        brick::numeric::Index2D keyPoint = keyPoints[ii];
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
  brick::computerVision::KeypointSelectorLepetitTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::KeypointSelectorLepetitTest currentTest;

}

#endif
