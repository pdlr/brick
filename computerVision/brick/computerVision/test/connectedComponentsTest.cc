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

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/connectedComponents.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/test/testFixture.hh>

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
      void testConnectedComponents();

    private:

    }; // class ConnectedComponentsTest


    /* ============== Member Function Definititions ============== */

    ConnectedComponentsTest::
    ConnectedComponentsTest()
      : TestFixture<ConnectedComponentsTest>("ConnectedComponentsTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testConnectedComponents);
    }


    void
    ConnectedComponentsTest::
    testConnectedComponents()
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
      Image<GRAY8> ccImage = connectedComponents<GRAY8>(binaryImage);
      for(size_t index0 = 0; index0 < inputImage.size(); ++index0) {
        BRICK_TEST_ASSERT(ccImage[index0] == inputImage[index0]);
      }
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
