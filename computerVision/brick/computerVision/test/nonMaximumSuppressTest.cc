/**
***************************************************************************
* @file nonMaximumSuppressTest.cc
*
* Source file defining tests for routines declared in
* brick/computerVision/nonMaximumSuppress.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/nonMaximumSuppress.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace computerVision {

    class NonMaximumSuppressTest
      : public brick::test::TestFixture<NonMaximumSuppressTest> {

    public:

      NonMaximumSuppressTest();
      ~NonMaximumSuppressTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testNonMaximumSuppress();

    private:

    }; // class NonMaximumSuppressTest


    /* ============== Member Function Definititions ============== */

    NonMaximumSuppressTest::
    NonMaximumSuppressTest()
      : brick::test::TestFixture<NonMaximumSuppressTest>("NonMaximumSuppressTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testNonMaximumSuppress);
    }


    void
    NonMaximumSuppressTest::
    testNonMaximumSuppress()
    {
      brick::numeric::Array2D<brick::common::Int32> inputArray(
        "[[ 0,  1,  0, -2, 5],"
        " [ 7,  5,  2,  3, 2],"
        " [ 8,  9,  8,  7, 6],"
        " [ 1, 10,  3,  4, 8],"
        " [ 2,  2,  5,  9, 2],"
        " [ 5,  4,  0,  6, 3],"
        " [ 1,  3,  5,  4, 3]]");


      brick::numeric::Array2D<brick::common::Float64> gradXArray(
        "[[ 1.0,  1.5,  1.7,  0.0, 1.0],"
        " [ 3.2,  0.0,  2.5,  6.1, 2.0],"
        " [-1.0,  9.0,  4.0,  1.2, 1.5],"
        " [-6.0,  1.5, -2.9, -6.0, 0.9],"
        " [-6.0, -1.9, -7.9, -5.0, 0.9],"
        " [-6.0, -0.7, -7.9,  4.0, 0.9],"
        " [   0,  1.2,  2.8,  3.2, 1.2]]");

      brick::numeric::Array2D<brick::common::Float64> gradYArray(
        "[[ 8.0,  0.0, -1.2,  -1.6,  1.2],"
        " [-1.8,  0.0,  0.9,   2.0, -1.0],"
        " [-1.0, -0.9, -1.4,  -0.7,  4.4],"
        " [ 2.6, -3.1, -9.3, -10.0, -2.1],"
        " [ 0.0, -2.0, -3.8,   1.0,  2.3],"
        " [ 0.0,  1.9, -3.8,   5.3,  2.3],"
        " [ 2.0,  2.5,  1.8,  -3.3,  1.0]]");

      brick::numeric::Array2D<brick::common::Float64> referenceArray(
        "[[  0,   0,  0,  0,  0],"
        " [  0,   0,  0,  3,  0],"
        " [  0,   9,  0,  7,  0],"
        " [  0,  10,  0,  0,  0],"
        " [  0,   2,  0,  9,  0],"
        " [  0,   4,  0,  6,  0],"
        " [  0,   0,  0,  0,  0]]");

      Image<GRAY_SIGNED32> inputImage = inputArray;
      Image<GRAY_FLOAT64> gradX = gradXArray;
      Image<GRAY_FLOAT64> gradY = gradYArray;
      Image<GRAY_SIGNED32> resultImage =
        nonMaximumSuppress<brick::common::Float64>(inputImage, gradX, gradY);

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
  brick::computerVision::NonMaximumSuppressTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::NonMaximumSuppressTest currentTest;

}

#endif
