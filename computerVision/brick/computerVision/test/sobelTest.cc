/**
***************************************************************************
* @file brick/computerVision/test/sobelTest.cc
*
* Source file defining tests for routines defined in
* brick/computerVision/sobel.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/sobel.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace computerVision {

    class SobelTest
      : public brick::test::TestFixture<SobelTest> {

    public:

      SobelTest();
      ~SobelTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testSobelX();
      void testSobelY();

    private:

    }; // class SobelTest


    /* ============== Member Function Definititions ============== */

    SobelTest::
    SobelTest()
      : brick::test::TestFixture<SobelTest>("SobelTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testSobelX);
      BRICK_TEST_REGISTER_MEMBER(testSobelY);
    }


    void
    SobelTest::
    testSobelX()
    {
      brick::numeric::Array2D<brick::common::Int32> inputArray(
        "[[ 0,  1,  0, -2],"
        " [ 7,  5,  2,  3],"
        " [10,  9,  8,  7],"
        " [ 3,  2,  3,  4],"
        " [ 2,  2,  5,  9]]");

      brick::numeric::Array2D<brick::common::Int32> referenceArray(
        "[[  8,   0, -12, -16],"
        " [ -8, -12,  -9,  -2],"
        " [-10,  -9,  -4,   0],"
        " [ -6,   1,   9,  10],"
        " [  0,  12,  28,  32]]");

      Image<GRAY_SIGNED32> inputImage = inputArray;
      Image<GRAY_SIGNED32> resultImage = applySobelX(inputImage);

      BRICK_TEST_ASSERT(resultImage.rows() == referenceArray.rows());
      BRICK_TEST_ASSERT(resultImage.columns() == referenceArray.columns());
      for(size_t index0 = 0; index0 < resultImage.size(); ++index0) {
        BRICK_TEST_ASSERT(resultImage[index0] == referenceArray[index0]);
      }
    }


    void
    SobelTest::
    testSobelY()
    {
      brick::numeric::Array2D<brick::common::Int32> inputArray(
        "[[ 0,  1,  0, -2],"
        " [ 7,  5,  2,  3],"
        " [10,  9,  8,  7],"
        " [ 3,  2,  3,  4],"
        " [ 2,  2,  5,  9]]");

      brick::numeric::Array2D<brick::common::Int32> referenceArray(
        "[[ 56, 34,  26,  40],"
        " [ 40, 34,  33,  36],"
        " [-16, -9,   0,   4],"
        " [-32,-25, -11,   8],"
        " [ -8,  2,  18,  40]]");

      Image<GRAY_SIGNED32> inputImage = inputArray;
      Image<GRAY_SIGNED32> resultImage = applySobelY(inputImage);

      BRICK_TEST_ASSERT(resultImage.rows() == referenceArray.rows());
      BRICK_TEST_ASSERT(resultImage.columns() == referenceArray.columns());
      for(size_t index0 = 0; index0 < resultImage.size(); ++index0) {
        BRICK_TEST_ASSERT(resultImage[index0] == referenceArray[index0]);
      }
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::SobelTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::SobelTest currentTest;

}

#endif
