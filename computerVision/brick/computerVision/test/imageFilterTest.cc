/**
***************************************************************************
* @file filterTest.cpp
*
* Source file defining tests for filter().
*
* Copyright (C) 2006 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/image.hh>
#include <brick/computerVision/imageFilter.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/kernel.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace computerVision {
    
    class FilterTest : public brick::test::TestFixture<FilterTest> {

    public:

      FilterTest();
      ~FilterTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testFilter2D_nonSeparable_i();
      void testFilter2D_nonSeparable();
      void testFilter2D_separable_i();
      void testFilter2D_separable();

    private:

      numeric::Array2D<common::UnsignedInt8>
      localFilter2D(const numeric::Array2D<common::Float64>& kernel,
                    const numeric::Array2D<common::UnsignedInt8>& inputImage);

    
      numeric::Array2D<common::Float64>
      localFilter2D(const numeric::Array2D<common::Float64>& kernel,
                    const numeric::Array2D<common::Float64>& inputImage);

      numeric::Array2D<common::UnsignedInt64>
      localFilter2D(const numeric::Array2D<common::UnsignedInt64>& kernel,
                    const numeric::Array2D<common::UnsignedInt8>& inputImage);

    }; // class FilterTest


    /* ============== Member Function Definititions ============== */

    FilterTest::
    FilterTest()
      : brick::test::TestFixture<FilterTest>("FilterTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testFilter2D_nonSeparable_i);
      BRICK_TEST_REGISTER_MEMBER(testFilter2D_nonSeparable);
      BRICK_TEST_REGISTER_MEMBER(testFilter2D_separable_i);
      BRICK_TEST_REGISTER_MEMBER(testFilter2D_separable);
    }


    void
    FilterTest::
    testFilter2D_nonSeparable_i()
    {
      Image<GRAY_FLOAT64> inputImage0(
	numeric::Array2D<common::Float64>("[[1.0, 2.0, 3.0, 4.0, 5.0],"
			 " [6.0, 7.0, 8.0, 9.0, 10.0],"
			 " [9.0, 8.0, 7.0, 6.0, 5.0],"
			 " [4.0, 3.0, 2.0, 1.0, 0.0]]"));
      numeric::Array2D<double> kernelData("[[1.0, 2.0, 1.0],"
                                 " [2.0, 4.0, 2.0],"
                                 " [3.0, 5.0, 4.0]]");
      // kernelData /= numeric::sum<double>(ravel(kernelData));
      Kernel<double> kernel0(kernelData);
      Image<GRAY_FLOAT64> resultImage =
	brick::computerVision::filter2D<
	GRAY_FLOAT64, GRAY_FLOAT64, GRAY_FLOAT64, common::Float64>(
	  kernel0, inputImage0, static_cast<common::Float64>(0));
      Image<GRAY_FLOAT64> referenceImage = localFilter2D(
        kernel0.getArray2D(), inputImage0);

      BRICK_TEST_ASSERT(resultImage.rows() == referenceImage.rows());
      BRICK_TEST_ASSERT(resultImage.columns() == referenceImage.columns());
      double tolerance = 1.0E-12;
      for(size_t index0 = 0; index0 < resultImage.size(); ++index0) {
        BRICK_TEST_ASSERT(
          approximatelyEqual(
            resultImage[index0], referenceImage[index0], tolerance));
      }
    }


    void
    FilterTest::
    testFilter2D_nonSeparable()
    {
      Image<GRAY8> inputImage0 = readPGM8(getTestImageFileNamePGM0());
      numeric::Array2D<double> kernelData("[[1.0, 2.0, 1.0],"
                                 " [2.0, 4.0, 2.0],"
                                 " [3.0, 5.0, 4.0],"
                                 " [2.0, 5.0, 1.0],"
                                 " [0.0, 2.0, 1.0]]");
      kernelData /= numeric::sum<double>(ravel(kernelData));
      Kernel<double> kernel0(kernelData);
      // Note(xxx): Need to make this not require explicit template
      // parameters.
      Image<GRAY8> resultImage =
        brick::computerVision::filter2D<GRAY8, GRAY_FLOAT64, GRAY8, double>(
          kernel0, inputImage0, static_cast<common::UnsignedInt8>(0));
      Image<GRAY8> referenceImage = localFilter2D(
        kernel0.getArray2D(), inputImage0);

      BRICK_TEST_ASSERT(resultImage.rows() == referenceImage.rows());
      BRICK_TEST_ASSERT(resultImage.columns() == referenceImage.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage.begin(), resultImage.end(),
                                 referenceImage.begin()));
    }


    void
    FilterTest::
    testFilter2D_separable_i()
    {
      Image<GRAY_FLOAT64> inputImage0(
	numeric::Array2D<common::Float64>("[[1.0, 2.0, 3.0, 4.0, 5.0],"
			 " [6.0, 7.0, 8.0, 9.0, 10.0],"
			 " [9.0, 8.0, 7.0, 6.0, 5.0],"
			 " [4.0, 3.0, 2.0, 1.0, 0.0]]"));

      numeric::Array1D<double> kernelRow("[1.0, 3.0, 2.0]");
      numeric::Array1D<double> kernelColumn("[2.0, 0.0, 1.0]");
      Kernel<double> kernel0(kernelRow, kernelColumn);

      numeric::Array2D<double> referenceKernelArray(
        kernelColumn.size(), kernelRow.size());
      for(size_t rowIndex = 0; rowIndex < referenceKernelArray.rows();
          ++rowIndex) {
        for(size_t columnIndex = 0;
	    columnIndex < referenceKernelArray.columns(); ++columnIndex) {
          referenceKernelArray(rowIndex, columnIndex) =
            kernelColumn[rowIndex] * kernelRow[columnIndex];
        }
      }


      Image<GRAY_FLOAT64> resultImage =
	brick::computerVision::filter2D<
	GRAY_FLOAT64, GRAY_FLOAT64, GRAY_FLOAT64, common::Float64>(
	  kernel0, inputImage0, static_cast<common::Float64>(0));
      Image<GRAY_FLOAT64> referenceImage = localFilter2D(
        kernel0.getArray2D(), inputImage0);

      BRICK_TEST_ASSERT(resultImage.rows() == referenceImage.rows());
      BRICK_TEST_ASSERT(resultImage.columns() == referenceImage.columns());
      double tolerance = 1.0E-12;
      for(size_t index0 = 0; index0 < resultImage.size(); ++index0) {
        BRICK_TEST_ASSERT(
          approximatelyEqual(
            resultImage[index0], referenceImage[index0], tolerance));
      }
    }


    void
    FilterTest::
    testFilter2D_separable()
    {
      Image<GRAY8> inputImage = readPGM8(getTestImageFileNamePGM0());
      Image<GRAY_FLOAT64> inputImage0(inputImage.rows(), inputImage.columns());
      inputImage0.copy(inputImage);
    
      numeric::Array1D<double> kernelRow("[1.0, 3.0, 4.0]");
      numeric::Array1D<double> kernelColumn("[3.0, 0.0, 2.0, 4.0, 1.0]");
      kernelRow /= numeric::sum<double>(kernelRow);
      kernelColumn /= numeric::sum<double>(kernelColumn);
      Kernel<double> kernel0(kernelRow, kernelColumn);

      numeric::Array2D<double> referenceKernelArray(
        kernelColumn.size(), kernelRow.size());
      for(size_t rowIndex = 0; rowIndex < referenceKernelArray.rows();
          ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < referenceKernelArray.columns();
            ++columnIndex) {
          referenceKernelArray(rowIndex, columnIndex) =
            kernelColumn[rowIndex] * kernelRow[columnIndex];
        }
      }

      // Note(xxx): Need to make this not require explicit template
      // parameters.
      Image<GRAY_FLOAT64> resultImage =
        computerVision::filter2D<
	GRAY_FLOAT64, GRAY_FLOAT64, GRAY_FLOAT64, double>(
          kernel0, inputImage0, static_cast<common::Float64>(0));
      Image<GRAY_FLOAT64> referenceImage = localFilter2D(
        referenceKernelArray, inputImage0);

      double tolerance = 1.0E-12;
      for(size_t index0 = 0; index0 < resultImage.size(); ++index0) {
        BRICK_TEST_ASSERT(
          approximatelyEqual(
            resultImage[index0], referenceImage[index0], tolerance));
      }
    }


    numeric::Array2D<common::UnsignedInt8>
    FilterTest::
    localFilter2D(const numeric::Array2D<common::Float64>& kernel,
                  const numeric::Array2D<common::UnsignedInt8>& inputImage)
    {
      size_t halfKernelRows = kernel.rows() / 2;
      size_t halfKernelColumns = kernel.columns() / 2;
      numeric::Array2D<common::UnsignedInt8> resultImage(inputImage.rows(), inputImage.columns());
      for(size_t imageRow = 0; imageRow < inputImage.rows(); ++imageRow) {
        for(size_t imageColumn = 0; imageColumn < inputImage.columns();
            ++imageColumn) {
          if(imageRow < halfKernelRows
             || imageColumn < halfKernelColumns
             || imageRow >= inputImage.rows() - halfKernelRows
             || imageColumn >= inputImage.columns() - halfKernelColumns) {
            resultImage(imageRow, imageColumn) = static_cast<common::UnsignedInt8>(0);
            continue;
          }
          common::Float64 dotProduct = 0.0;
          for(size_t kernelRow = 0; kernelRow < kernel.rows(); ++kernelRow) {
            for(size_t kernelColumn = 0; kernelColumn < kernel.columns();
                ++kernelColumn) {
              dotProduct += (
                inputImage(imageRow - halfKernelRows + kernelRow,
                           imageColumn - halfKernelColumns + kernelColumn)
                * kernel(kernelRow, kernelColumn));
            }
          }
          resultImage(imageRow, imageColumn) =
//             static_cast<common::UnsignedInt8>(dotProduct + 0.5);
	    static_cast<common::UnsignedInt8>(dotProduct);
        }
      }
      return resultImage;
    }


    numeric::Array2D<common::Float64>
    FilterTest::
    localFilter2D(const numeric::Array2D<common::Float64>& kernel,
                  const numeric::Array2D<common::Float64>& inputImage)
    {
      size_t halfKernelRows = kernel.rows() / 2;
      size_t halfKernelColumns = kernel.columns() / 2;
      numeric::Array2D<common::Float64> resultImage(inputImage.rows(), inputImage.columns());
      for(size_t imageRow = 0; imageRow < inputImage.rows(); ++imageRow) {
        for(size_t imageColumn = 0; imageColumn < inputImage.columns();
            ++imageColumn) {
          if(imageRow < halfKernelRows
             || imageColumn < halfKernelColumns
             || imageRow >= inputImage.rows() - halfKernelRows
             || imageColumn >= inputImage.columns() - halfKernelColumns) {
            resultImage(imageRow, imageColumn) = static_cast<common::Float64>(0);
            continue;
          }
          common::Float64 dotProduct = 0.0;
          for(size_t kernelRow = 0; kernelRow < kernel.rows(); ++kernelRow) {
            for(size_t kernelColumn = 0; kernelColumn < kernel.columns();
                ++kernelColumn) {
              dotProduct += (
                inputImage(imageRow - halfKernelRows + kernelRow,
                           imageColumn - halfKernelColumns + kernelColumn)
                * kernel(kernelRow, kernelColumn));
            }
          }
          resultImage(imageRow, imageColumn) =
            static_cast<common::Float64>(dotProduct);
        }
      }
      return resultImage;
    }  


    numeric::Array2D<common::UnsignedInt64>
    FilterTest::
    localFilter2D(const numeric::Array2D<common::UnsignedInt64>& kernel,
                  const numeric::Array2D<common::UnsignedInt8>& inputImage)
    {
      size_t halfKernelRows = kernel.rows() / 2;
      size_t halfKernelColumns = kernel.columns() / 2;
      numeric::Array2D<common::UnsignedInt64> resultImage(inputImage.rows(), inputImage.columns());
      for(size_t imageRow = 0; imageRow < inputImage.rows(); ++imageRow) {
        for(size_t imageColumn = 0; imageColumn < inputImage.columns();
            ++imageColumn) {
          if(imageRow < halfKernelRows
             || imageColumn < halfKernelColumns
             || imageRow >= inputImage.rows() - halfKernelRows
             || imageColumn >= inputImage.columns() - halfKernelColumns) {
            resultImage(imageRow, imageColumn) = static_cast<common::UnsignedInt64>(0);
            continue;
          }
          common::UnsignedInt64 dotProduct = 0;
          for(size_t kernelRow = 0; kernelRow < kernel.rows(); ++kernelRow) {
            for(size_t kernelColumn = 0; kernelColumn < kernel.columns();
                ++kernelColumn) {
              dotProduct += (
                inputImage(imageRow - halfKernelRows + kernelRow,
                           imageColumn - halfKernelColumns + kernelColumn)
                * kernel(kernelRow, kernelColumn));
            }
          }
          resultImage(imageRow, imageColumn) = dotProduct;
        }
      }
      return resultImage;
    }  

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::FilterTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::FilterTest currentTest;

}

#endif
