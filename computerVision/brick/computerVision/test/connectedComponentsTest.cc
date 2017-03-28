/**
***************************************************************************
* @file brick/computerVision/test/connectedComponentsTest.cc
*
* Source file defining tests for connectedComponents().
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <set>

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/connectedComponents.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/test/testFixture.hh>

#include <brick/portability/timeUtilities.hh>

namespace brick {

  namespace computerVision {

    class ConnectedComponentsTest
      : public brick::test::TestFixture<ConnectedComponentsTest> {

    public:

      ConnectedComponentsTest();
      ~ConnectedComponentsTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testConnectedComponentsForegroundBackground();
      void testConnectedComponentsSameColor();
      void testConnectedComponentsTiming();

    private:

    }; // class ConnectedComponentsTest


    /* ============== Member Function Definititions ============== */

    ConnectedComponentsTest::
    ConnectedComponentsTest()
      : brick::test::TestFixture<ConnectedComponentsTest>("ConnectedComponentsTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testConnectedComponentsForegroundBackground);
      BRICK_TEST_REGISTER_MEMBER(testConnectedComponentsSameColor);
      // BRICK_TEST_REGISTER_MEMBER(testConnectedComponentsTiming);
    }


    void
    ConnectedComponentsTest::
    testConnectedComponentsForegroundBackground()
    {
      // Input image has labled connected components.
      Image<GRAY8> inputImage = readPGM8(getConnectedComponentsFileNamePGM0());

      // Generate an image where all foreground components are
      // nonzero, and all background components are zero.
      Image<GRAY8> binaryImage(inputImage.rows(), inputImage.columns());
      for(size_t index0 = 0; index0 < inputImage.size(); ++index0) {
        if(inputImage[index0] == 0) {
          binaryImage[index0] = brick::common::UnsignedInt8(0);
        } else {
          binaryImage[index0] = brick::common::UnsignedInt8(1);
        }
      }

      // Run connected components in binary mode.
      unsigned int numberOfComponents = 0;
      Image<GRAY8> ccImage = connectedComponents<GRAY8>(
        binaryImage, numberOfComponents);

      // Make sure we were able to recover the correct labels.
      // This test assumes (because of implementation knowledge) that
      // labels increase as objects are encountered in raster order.
      // Note that this assumption is not guaranteed by the published
      // interface.
      for(size_t index0 = 0; index0 < inputImage.size(); ++index0) {
        BRICK_TEST_ASSERT(ccImage[index0] == inputImage[index0]);
      }

      // Make sure we reported the right number of components.
      std::set<uint8_t> labelSet;
      for(auto ccIter = ccImage.begin(); ccIter != ccImage.end(); ++ccIter) {
        if(*ccIter != 0) {
          labelSet.insert(*ccIter);
        }
      }
      BRICK_TEST_ASSERT(numberOfComponents == labelSet.size());
    }

    void
    ConnectedComponentsTest::
    testConnectedComponentsSameColor()
    {
      // Input image has labled connected components, but the labels
      // are ground truth for connected components analysis in binary
      // mode.  No matter, we can threshold this to make a good input
      // for our non-binary connected components routine.
      Image<GRAY8> inputImage = readPGM8(
        getConnectedComponentsFileNamePGM0());
      Image<GRAY8> binaryImage(inputImage.rows(), inputImage.columns());
      for(size_t index0 = 0; index0 < inputImage.size(); ++index0) {
        if(inputImage[index0] == 0) {
          binaryImage[index0] = brick::common::UnsignedInt8(0);
        } else {
          binaryImage[index0] = brick::common::UnsignedInt8(1);
        }
      }
      
      // This is the same image, but with labels that reflect ground
      // truth for connected components run in non-binary mode.
      Image<GRAY8> groundTruthImage = readPGM8(
        getConnectedComponentsFileNamePGM1());

      // Run connected components in non-binary mode.
      ConnectedComponentsConfig config;
      config.mode = ConnectedComponentsConfig::SAME_COLOR;
      unsigned int numberOfComponents = 0;
      Image<GRAY8> ccImage = connectedComponents<GRAY8>(
        binaryImage, numberOfComponents, config);

      // Make sure we were able to recover the correct labels.
      // This test assumes (because of implementation knowledge) that
      // labels increase as objects are encountered in raster order.
      // Note that this assumption is not guaranteed by the published
      // interface.
      for(size_t index0 = 0; index0 < inputImage.size(); ++index0) {
        BRICK_TEST_ASSERT(ccImage[index0] == groundTruthImage[index0]);
      }

      // Make sure we reported the right number of components.
      std::set<uint8_t> labelSet;
      for(auto ccIter = ccImage.begin(); ccIter != ccImage.end(); ++ccIter) {
        labelSet.insert(*ccIter);
      }
      BRICK_TEST_ASSERT(numberOfComponents == labelSet.size());
    }

    void
    ConnectedComponentsTest::
    testConnectedComponentsTiming()
    {
      Image<GRAY8> inputImage = readPGM8(getConnectedComponentsFileNamePGM0());
      Image<GRAY8> binaryImage(inputImage.rows(), inputImage.columns());
       for(size_t index0 = 0; index0 < inputImage.size(); ++index0) {
        if(inputImage[index0] == 0) {
          binaryImage[index0] = brick::common::UnsignedInt8(0);
        } else {
          binaryImage[index0] = brick::common::UnsignedInt8(1);
        }
      }

      Image<GRAY8> ccImage;
      unsigned int numberOfComponents = 0;

      ConnectedComponentsConfig config0;
      config0.mode = ConnectedComponentsConfig::FOREGROUND_BACKGROUND;
      ConnectedComponentsConfig config1;
      config1.mode = ConnectedComponentsConfig::SAME_COLOR;

      double t0 = brick::portability::getCurrentTime();
      for(int ii = 0; ii < 100; ++ii) {
        ccImage = connectedComponents<GRAY8>(binaryImage, numberOfComponents,
                                             config0);
      }
      double t1 = brick::portability::getCurrentTime();
      for(int ii = 0; ii < 100; ++ii) {
        ccImage = connectedComponents<GRAY8>(binaryImage, numberOfComponents,
                                             config1);
      }
      double t2 = brick::portability::getCurrentTime();
      for(int ii = 0; ii < 100; ++ii) {
        ccImage = connectedComponents<GRAY8>(binaryImage, numberOfComponents,
                                             config0);
      }
      double t3 = brick::portability::getCurrentTime();
      for(int ii = 0; ii < 100; ++ii) {
        ccImage = connectedComponents<GRAY8>(binaryImage, numberOfComponents,
                                             config1);
      }
      double t4 = brick::portability::getCurrentTime();

      std::cout << "Config0: " << t1 - t0 << std::endl;
      std::cout << "Config1: " << t2 - t1 << std::endl;
      std::cout << "Config0: " << t3 - t2 << std::endl;
      std::cout << "Config1: " << t4 - t3 << std::endl;
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::ConnectedComponentsTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ConnectedComponentsTest currentTest;

}

#endif
