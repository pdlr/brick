/**
***************************************************************************
* @file thresholderSauvolaTest.cpp
*
* Source file defining tests for the ThresholderSauvola class.
*
* Copyright (C) 2013 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <iomanip>
#include <iostream>

#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/thresholderSauvola.hh>
#include <brick/computerVision/utilities.hh>
#include <brick/computerVision/test/testImages.hh>

#include <brick/numeric/subArray2D.hh>

#include <brick/test/testFixture.hh>

#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace computerVision {


    class ThresholderSauvolaTest
      : public brick::test::TestFixture<ThresholderSauvolaTest> {

    public:

      ThresholderSauvolaTest();
      ~ThresholderSauvolaTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testThresholderSauvola();
      void testExecutionTime();

    private:

      double m_defaultTolerance;
      uint32_t m_kernelSize;

    }; // class ThresholderSauvolaTest


    /* ============== Member Function Definititions ============== */

    ThresholderSauvolaTest::
    ThresholderSauvolaTest()
      : brick::test::TestFixture<ThresholderSauvolaTest>(
          "ThresholderSauvolaTest"),
        m_defaultTolerance(1.0E-8),
        m_kernelSize(64)
    {
      BRICK_TEST_REGISTER_MEMBER(testThresholderSauvola);
      // BRICK_TEST_REGISTER_MEMBER(testExecutionTime);
    }


    void
    ThresholderSauvolaTest::
    testThresholderSauvola()
    {
      // The input file ought to be text, but we'll use a natural
      // scene to get started.
      Image<GRAY8> inputImage = readPGM8(getBullseyeFileNamePGM0());

      // For now, just threshold the image and see how it looks.
      ThresholderSauvola<GRAY8> thresholder(m_kernelSize, 0.5);
      thresholder.setImage(inputImage);
      Image<GRAY8> outputImage = thresholder.computeBinaryImage();

      BRICK_TEST_ASSERT(outputImage.rows() == inputImage.rows());
      BRICK_TEST_ASSERT(outputImage.columns() == inputImage.columns());

      // writePGM8("foo.pgm", outputImage);
    }


    void
    ThresholderSauvolaTest::
    testExecutionTime()
    {
      // Execution time doesn't depend on image content, so we'll use
      // whatever image is handy.
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

      // Run the thresholder.
      ThresholderSauvola<GRAY8> thresholder(m_kernelSize, 0.5);

      std::cout << "Start!" << std::endl;

      double t0 = brick::utilities::getCurrentTime();
      thresholder.setImage(bigImage);
      Image<GRAY8> outputImage = thresholder.computeBinaryImage();
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
  brick::computerVision::ThresholderSauvolaTest currentTest;
  // currentTest.testThresholderSauvola();
  bool result = currentTest.run();

  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ThresholderSauvolaTest currentTest;

}

#endif
