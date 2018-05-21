/**
***************************************************************************
* @file brick/numeric/test/stencil2DTest.cpp
* 
* Source file defining Stencil2DTest class.
*
* Copyright (C) 2004-2005 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/stencil2D.hh>
#include <brick/portability/timeUtilities.hh>
#include <brick/test/testFixture.hh>

// using namespace brick::numeric;
using namespace brick;
using brick::portability::getCurrentTime;

namespace brick {

  namespace numeric {
    
    class Stencil2DTest
      : public brick::test::TestFixture<Stencil2DTest> {

    public:

      Stencil2DTest();
      ~Stencil2DTest() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testSpeed();
    
    private:

      void
      doAllOutConvolution(const Array2D<int>& inputImage,
                          const Array2D<int>& kernel,
                          Array2D<int>& resultImage);

    
      void
      doArrayDerefConvolution(const Array2D<int>& inputImage,
                              const Array2D<int>& kernel,
                              Array2D<int>& resultImage);

    

      void
      doConstArrayDerefConvolution(const Array2D<int>& inputImage,
                                   const Array2D<int>& kernel,
                                   Array2D<int>& resultImage);

    
      void
      doConstVariableDerefConvolution(const Array2D<int>& inputImage,
                                      const Array2D<int>& kernel,
                                      Array2D<int>& resultImage);

    
      void
      doDlrArrayDerefConvolution(const Array2D<int>& inputImage,
                                 const Array2D<int>& kernel,
                                 Array2D<int>& resultImage);

    

      void
      doHalfIteratorConvolution(const Array2D<int>& inputImage,
                                const Array2D<int>& kernel,
                                Array2D<int>& resultImage);
    

    
      void
      doHeapArrayDerefConvolution(const Array2D<int>& inputImage,
                                  const Array2D<int>& kernel,
                                  Array2D<int>& resultImage);

    
      void
      doIteratorConvolution(const Array2D<int>& inputImage,
                            const Array2D<int>& kernel,
                            Array2D<int>& resultImage);

    
      void
      doNaiveConvolution(const Array2D<int>& inputImage,
                         const Array2D<int>& kernel,
                         Array2D<int>& resultImage);
    

      void
      doRunTimeStencilConvolution(const Array2D<int>& inputImage,
                                  const Array2D<int>& kernel,
                                  Array2D<int>& resultImage);

    
      void
      doSimplifiedStencilConvolution(const Array2D<int>& inputImage,
                                     const Array2D<int>& kernel,
                                     Array2D<int>& resultImage);


      void
      doSingleIndexConvolution(const Array2D<int>& inputImage,
                               const Array2D<int>& kernel,
                               Array2D<int>& resultImage);
    

      void
      doStructDerefConvolution(const Array2D<int>& inputImage,
                               const Array2D<int>& kernel,
                               Array2D<int>& resultImage);
    
      
      void
      doStructArrayDerefConvolution(const Array2D<int>& inputImage,
                                    const Array2D<int>& kernel,
                                    Array2D<int>& resultImage);
    
      
      void
      doVariableDerefConvolution(const Array2D<int>& inputImage,
                                 const Array2D<int>& kernel,
                                 Array2D<int>& resultImage);

    
    }; // class Stencil2DTest


    /* ============== Member Function Definititions ============== */

    Stencil2DTest::
    Stencil2DTest()
      : brick::test::TestFixture<Stencil2DTest>("Stencil2DTest")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testSpeed);
    }


    void
    Stencil2DTest::
    testSpeed()
    {
      // Do a 3x3 convolution in several different ways to test speed.
      const size_t imageSize = 512;
      const size_t kernelSize = 3;
      const int imageValue = 3;
      const int kernelValue = 4;
    
      Array2D<int> image0(imageSize, imageSize);
      Array2D<int> kernel(kernelSize, kernelSize);

      for(size_t rowIndex = 0; rowIndex < image0.rows(); ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < image0.columns();
            ++columnIndex) {
          image0(rowIndex, columnIndex) = static_cast<int>((rowIndex + 2 * columnIndex) % 124);
        }
      }

      int count = 0;
      for(size_t rowIndex = 0; rowIndex < kernel.rows(); ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < kernel.columns();
            ++columnIndex) {
          image0 = static_cast<int>(kernel.size()) - count;
        }
      }

      Array2D<int> resultImage0(imageSize, imageSize);
      Array2D<int> resultImage1(imageSize, imageSize);

      image0 = imageValue;
      kernel = kernelValue;
    
      resultImage0 = 0;
      double t0 = getCurrentTime();
      this->doNaiveConvolution(image0, kernel, resultImage0);
      double t1 = getCurrentTime();
      std::cout << "Reference convolution: " << t1 - t0 << std::endl;

      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doNaiveConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Naive convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doSingleIndexConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Single index convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doRunTimeStencilConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Run-time stencil convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doAllOutConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "All-out convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doConstVariableDerefConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Const variable deref convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doVariableDerefConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Variable deref convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doConstArrayDerefConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Const array deref convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doArrayDerefConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Array deref convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doDlrArrayDerefConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Dlr array deref convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doHeapArrayDerefConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Heap array deref convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doIteratorConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Iterator convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doHalfIteratorConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Half Iterator convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doStructDerefConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Struct deref convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doStructArrayDerefConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Struct array deref convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));


      resultImage1 = 0;
      t0 = getCurrentTime();
      this->doSimplifiedStencilConvolution(image0, kernel, resultImage1);
      t1 = getCurrentTime();
      std::cout << "Simplified stencil convolution: " << t1 - t0 << std::endl;
      BRICK_TEST_ASSERT(resultImage1.rows() == resultImage0.rows());
      BRICK_TEST_ASSERT(resultImage1.columns() == resultImage0.columns());
      BRICK_TEST_ASSERT(std::equal(resultImage1.begin(), resultImage1.end(),
                                   resultImage0.begin()));



    }


    void
    Stencil2DTest::
    doAllOutConvolution(const Array2D<int>& inputImage,
                        const Array2D<int>& kernel,
                        Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();
    
      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          dotProduct += *(kernelPtr + 0) * *(imagePtr + 0);
          dotProduct += *(kernelPtr + 1) * *(imagePtr + 1);
          dotProduct += *(kernelPtr + 2) * *(imagePtr + 2);
          dotProduct += *(kernelPtr + 3) * *(imagePtr + 512);
          dotProduct += *(kernelPtr + 4) * *(imagePtr + 513);
          dotProduct += *(kernelPtr + 5) * *(imagePtr + 514);
          dotProduct += *(kernelPtr + 6) * *(imagePtr + 1024);
          dotProduct += *(kernelPtr + 7) * *(imagePtr + 1025);
          dotProduct += *(kernelPtr + 8) * *(imagePtr + 1026);
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }
    }


    void
    Stencil2DTest::
    doArrayDerefConvolution(const Array2D<int>& inputImage,
                            const Array2D<int>& kernel,
                            Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();

      int ii[9] = {0,  1,  2,  512,  513,  514,  1024,  1025,  1026};
      int kk[9] = {0,  1,  2,  3,  4,  5,  6,  7,  8};

      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          dotProduct += *(kernelPtr + kk[0]) * *(imagePtr + ii[0]);
          dotProduct += *(kernelPtr + kk[1]) * *(imagePtr + ii[1]);
          dotProduct += *(kernelPtr + kk[2]) * *(imagePtr + ii[2]);
          dotProduct += *(kernelPtr + kk[3]) * *(imagePtr + ii[3]);
          dotProduct += *(kernelPtr + kk[4]) * *(imagePtr + ii[4]);
          dotProduct += *(kernelPtr + kk[5]) * *(imagePtr + ii[5]);
          dotProduct += *(kernelPtr + kk[6]) * *(imagePtr + ii[6]);
          dotProduct += *(kernelPtr + kk[7]) * *(imagePtr + ii[7]);
          dotProduct += *(kernelPtr + kk[8]) * *(imagePtr + ii[8]);
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }
    }


    void
    Stencil2DTest::
    doConstArrayDerefConvolution(const Array2D<int>& inputImage,
                                 const Array2D<int>& kernel,
                                 Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();

      const int ii[9] = {0,  1,  2,  512,  513,  514,  1024,  1025,  1026};
      const int kk[9] = {0,  1,  2,  3,  4,  5,  6,  7,  8};

      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          dotProduct += *(kernelPtr + kk[0]) * *(imagePtr + ii[0]);
          dotProduct += *(kernelPtr + kk[1]) * *(imagePtr + ii[1]);
          dotProduct += *(kernelPtr + kk[2]) * *(imagePtr + ii[2]);
          dotProduct += *(kernelPtr + kk[3]) * *(imagePtr + ii[3]);
          dotProduct += *(kernelPtr + kk[4]) * *(imagePtr + ii[4]);
          dotProduct += *(kernelPtr + kk[5]) * *(imagePtr + ii[5]);
          dotProduct += *(kernelPtr + kk[6]) * *(imagePtr + ii[6]);
          dotProduct += *(kernelPtr + kk[7]) * *(imagePtr + ii[7]);
          dotProduct += *(kernelPtr + kk[8]) * *(imagePtr + ii[8]);
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }
    }


    void
    Stencil2DTest::
    doConstVariableDerefConvolution(const Array2D<int>& inputImage,
                                    const Array2D<int>& kernel,
                                    Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();

      const int i0 = 0;
      const int i1 = 1;
      const int i2 = 2;
      const int i3 = 512;
      const int i4 = 513;
      const int i5 = 514;
      const int i6 = 1024;
      const int i7 = 1025;
      const int i8 = 1026;

      const int k0 = 0;
      const int k1 = 1;
      const int k2 = 2;
      const int k3 = 3;
      const int k4 = 4;
      const int k5 = 5;
      const int k6 = 6;
      const int k7 = 7;
      const int k8 = 8;

      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          dotProduct += *(kernelPtr + k0) * *(imagePtr + i0);
          dotProduct += *(kernelPtr + k1) * *(imagePtr + i1);
          dotProduct += *(kernelPtr + k2) * *(imagePtr + i2);
          dotProduct += *(kernelPtr + k3) * *(imagePtr + i3);
          dotProduct += *(kernelPtr + k4) * *(imagePtr + i4);
          dotProduct += *(kernelPtr + k5) * *(imagePtr + i5);
          dotProduct += *(kernelPtr + k6) * *(imagePtr + i6);
          dotProduct += *(kernelPtr + k7) * *(imagePtr + i7);
          dotProduct += *(kernelPtr + k8) * *(imagePtr + i8);
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }
    }


    void
    Stencil2DTest::
    doDlrArrayDerefConvolution(const Array2D<int>& inputImage,
                               const Array2D<int>& kernel,
                               Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();

      Array1D<int> iiii(9);
      Array1D<int> kkkk(9);
      iiii[0] = 0;
      iiii[1] = 1;
      iiii[2] = 2;
      iiii[3] = 512;
      iiii[4] = 513;
      iiii[5] = 514;
      iiii[6] = 1024;
      iiii[7] = 1025;
      iiii[8] = 1026;

      kkkk[0] = 0;
      kkkk[1] = 1;
      kkkk[2] = 2;
      kkkk[3] = 3;
      kkkk[4] = 4;
      kkkk[5] = 5;
      kkkk[6] = 6;
      kkkk[7] = 7;
      kkkk[8] = 8;

      int* ii = iiii.data();
      int* kk = kkkk.data();
    
      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          dotProduct += *(kernelPtr + kk[0]) * *(imagePtr + ii[0]);
          dotProduct += *(kernelPtr + kk[1]) * *(imagePtr + ii[1]);
          dotProduct += *(kernelPtr + kk[2]) * *(imagePtr + ii[2]);
          dotProduct += *(kernelPtr + kk[3]) * *(imagePtr + ii[3]);
          dotProduct += *(kernelPtr + kk[4]) * *(imagePtr + ii[4]);
          dotProduct += *(kernelPtr + kk[5]) * *(imagePtr + ii[5]);
          dotProduct += *(kernelPtr + kk[6]) * *(imagePtr + ii[6]);
          dotProduct += *(kernelPtr + kk[7]) * *(imagePtr + ii[7]);
          dotProduct += *(kernelPtr + kk[8]) * *(imagePtr + ii[8]);
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }
    }


    void
    Stencil2DTest::
    doHalfIteratorConvolution(const Array2D<int>& inputImage,
                              const Array2D<int>& kernel,
                              Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      Stencil2D<const int, 9> imageStencil(kernelSize, kernelSize);
      imageStencil.setTarget(inputImage);

      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        imageStencil.goTo(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          Array2D<int>::const_iterator kernelIter = kernel.begin();
          Array2D<int>::const_iterator endIter = kernel.end();
          StencilIterator<const int, 9> imageIter = imageStencil.begin();
          while(kernelIter != endIter) {
            dotProduct += *kernelIter * *imageIter;
            ++kernelIter;
            ++imageIter;
          }
          resultImage(resultIndex) = dotProduct;
          ++resultIndex;
          imageStencil.advance();
        }
      }
    }


    void
    Stencil2DTest::
    doHeapArrayDerefConvolution(const Array2D<int>& inputImage,
                                const Array2D<int>& kernel,
                                Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();

      int* ii = new int[9];
      int* kk = new int[9];
      ii[0] = 0;
      ii[1] = 1;
      ii[2] = 2;
      ii[3] = 512;
      ii[4] = 513;
      ii[5] = 514;
      ii[6] = 1024;
      ii[7] = 1025;
      ii[8] = 1026;

      kk[0] = 0;
      kk[1] = 1;
      kk[2] = 2;
      kk[3] = 3;
      kk[4] = 4;
      kk[5] = 5;
      kk[6] = 6;
      kk[7] = 7;
      kk[8] = 8;

      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          for(size_t index0 = 0; index0 < 9; ++index0) {
            dotProduct += *(kernelPtr + kk[index0]) * *(imagePtr + ii[index0]);
          }
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }

      delete[] ii;
      delete[] kk;
    }


    void
    Stencil2DTest::
    doIteratorConvolution(const Array2D<int>& inputImage,
                          const Array2D<int>& kernel,
                          Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      Stencil2D<const int, 9> imageStencil(kernelSize, kernelSize);
      imageStencil.setTarget(inputImage);

      Stencil2D<const int, 9> kernelStencil(kernelSize, kernelSize);
      kernelStencil.setTarget(kernel);

      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        imageStencil.goTo(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          StencilIterator<const int, 9> kernelIter = kernelStencil.begin();
          StencilIterator<const int, 9> endIter = kernelStencil.end();
          StencilIterator<const int, 9> imageIter = imageStencil.begin();
          while(kernelIter != endIter) {
            dotProduct += *kernelIter * *imageIter;
            ++kernelIter;
            ++imageIter;
          }
          resultImage(resultIndex) = dotProduct;
          ++resultIndex;
          imageStencil.advance();
        }
      }
    }


    void
    Stencil2DTest::
    doNaiveConvolution(const Array2D<int>& inputImage,
                       const Array2D<int>& kernel,
                       Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t kernelRadius = kernelSize / 2;
      const size_t imageSize = inputImage.rows();

      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;
    
      for(size_t row = startRow; row < stopRow; ++row) {
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          for(size_t kernelRow = 0; kernelRow < kernelSize;
              ++kernelRow) {
            for(size_t kernelColumn = 0; kernelColumn < kernelSize;
                ++kernelColumn) {
              dotProduct +=
                inputImage(row, column) * kernel(kernelRow, kernelColumn);
            }
          }
          resultImage(row + kernelRadius, column + kernelRadius) = dotProduct;
        }
      }
    }
  

    void
    Stencil2DTest::
    doRunTimeStencilConvolution(const Array2D<int>& inputImage,
                                const Array2D<int>& kernel,
                                Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      Stencil2D<const int, 9> imageStencil(kernelSize, kernelSize);
      imageStencil.setTarget(inputImage);

      Stencil2D<const int, 9> kernelStencil(kernelSize, kernelSize);
      kernelStencil.setTarget(kernel);

      size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        imageStencil.goTo(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          for(size_t stencilIndex = 0; stencilIndex < stencilSize;
              ++stencilIndex) {
            dotProduct += (kernelStencil.getValue(stencilIndex)
                           * imageStencil.getValue(stencilIndex));
          }
          resultImage(resultIndex) = dotProduct;
          ++resultIndex;
          imageStencil.advance();
        }
      }
    }


    void
    Stencil2DTest::
    doSimplifiedStencilConvolution(const Array2D<int>& inputImage,
                                   const Array2D<int>& kernel,
                                   Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();

      Stencil2D<const int, 9> ii(kernelSize, kernelSize);
      Stencil2D<const int, 9> kk(kernelSize, kernelSize);
      ii.setTarget(inputImage);
      kk.setTarget(kernel);
    
      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          for(size_t index0 = 0; index0 < 9; ++index0) {
            dotProduct +=
              *(kernelPtr + kk.m_offsetArray[index0])
              * *(imagePtr + ii.m_offsetArray[index0]);
          }
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }
    }


    void
    Stencil2DTest::
    doSingleIndexConvolution(const Array2D<int>& inputImage,
                             const Array2D<int>& kernel,
                             Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;
    
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = row * imageSize + startColumn;
        for(size_t column = startColumn; column < stopColumn; ++column) {
          size_t imageIndex = resultIndex;
          size_t kernelIndex = 0;
          int dotProduct = 0;
          for(size_t kernelRow = 0; kernelRow < kernelSize;
              ++kernelRow) {
            for(size_t kernelColumn = 0; kernelColumn < kernelSize;
                ++kernelColumn) {
              dotProduct += inputImage(imageIndex) * kernel(kernelIndex);
              ++imageIndex;
              ++kernelIndex;
            }
            imageIndex += imageSize - kernelSize;
          }
          resultImage(resultIndex + imageSize + 1) = dotProduct;
          ++resultIndex;
        }
      }
    }


    struct StencilStruct {
      int c0;
      int c1;
      int c2;
      int c3;
      int c4;
      int c5;
      int c6;
      int c7;
      int c8;
      int cc[9];
    };

  
    void
    Stencil2DTest::
    doStructDerefConvolution(const Array2D<int>& inputImage,
                             const Array2D<int>& kernel,
                             Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();

      StencilStruct ii;
      StencilStruct kk;
      ii.c0 = 0;
      ii.c1 = 1;
      ii.c2 = 2;
      ii.c3 = 512;
      ii.c4 = 513;
      ii.c5 = 514;
      ii.c6 = 1024;
      ii.c7 = 1025;
      ii.c8 = 1026;

      kk.c0 = 0;
      kk.c1 = 1;
      kk.c2 = 2;
      kk.c3 = 3;
      kk.c4 = 4;
      kk.c5 = 5;
      kk.c6 = 6;
      kk.c7 = 7;
      kk.c8 = 8;

      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          dotProduct += *(kernelPtr + kk.c0) * *(imagePtr + ii.c0);
          dotProduct += *(kernelPtr + kk.c1) * *(imagePtr + ii.c1);
          dotProduct += *(kernelPtr + kk.c2) * *(imagePtr + ii.c2);
          dotProduct += *(kernelPtr + kk.c3) * *(imagePtr + ii.c3);
          dotProduct += *(kernelPtr + kk.c4) * *(imagePtr + ii.c4);
          dotProduct += *(kernelPtr + kk.c5) * *(imagePtr + ii.c5);
          dotProduct += *(kernelPtr + kk.c6) * *(imagePtr + ii.c6);
          dotProduct += *(kernelPtr + kk.c7) * *(imagePtr + ii.c7);
          dotProduct += *(kernelPtr + kk.c8) * *(imagePtr + ii.c8);
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }
    }


    void
    Stencil2DTest::
    doStructArrayDerefConvolution(const Array2D<int>& inputImage,
                                  const Array2D<int>& kernel,
                                  Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();

      StencilStruct ii;
      StencilStruct kk;
      ii.cc[0] = 0;
      ii.cc[1] = 1;
      ii.cc[2] = 2;
      ii.cc[3] = 512;
      ii.cc[4] = 513;
      ii.cc[5] = 514;
      ii.cc[6] = 1024;
      ii.cc[7] = 1025;
      ii.cc[8] = 1026;

      kk.cc[0] = 0;
      kk.cc[1] = 1;
      kk.cc[2] = 2;
      kk.cc[3] = 3;
      kk.cc[4] = 4;
      kk.cc[5] = 5;
      kk.cc[6] = 6;
      kk.cc[7] = 7;
      kk.cc[8] = 8;

      size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          for(size_t index0 = 0; index0 < stencilSize; ++index0) {
            dotProduct +=
              *(kernelPtr + kk.cc[index0]) * *(imagePtr + ii.cc[index0]);
          }
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }
    }


    void
    Stencil2DTest::
    doVariableDerefConvolution(const Array2D<int>& inputImage,
                               const Array2D<int>& kernel,
                               Array2D<int>& resultImage)
    {
      const size_t kernelSize = kernel.rows();
      const size_t imageSize = inputImage.rows();

      BRICK_TEST_ASSERT(kernelSize == 3);
      
      const size_t startRow = 0;
      const size_t stopRow = imageSize - kernelSize;
      const size_t startColumn = 0;
      const size_t stopColumn = imageSize - kernelSize;

      const int* kernelPtr = kernel.data();

      int i0 = 0;
      int i1 = 1;
      int i2 = 2;
      int i3 = 512;
      int i4 = 513;
      int i5 = 514;
      int i6 = 1024;
      int i7 = 1025;
      int i8 = 1026;

      int k0 = 0;
      int k1 = 1;
      int k2 = 2;
      int k3 = 3;
      int k4 = 4;
      int k5 = 5;
      int k6 = 6;
      int k7 = 7;
      int k8 = 8;

      // size_t stencilSize = kernelSize * kernelSize;
      for(size_t row = startRow; row < stopRow; ++row) {
        size_t resultIndex = (row + 1) * imageSize + startColumn + 1;
        const int* imagePtr = inputImage.data(row, startColumn);
        for(size_t column = startColumn; column < stopColumn; ++column) {
          int dotProduct = 0;
          dotProduct += *(kernelPtr + k0) * *(imagePtr + i0);
          dotProduct += *(kernelPtr + k1) * *(imagePtr + i1);
          dotProduct += *(kernelPtr + k2) * *(imagePtr + i2);
          dotProduct += *(kernelPtr + k3) * *(imagePtr + i3);
          dotProduct += *(kernelPtr + k4) * *(imagePtr + i4);
          dotProduct += *(kernelPtr + k5) * *(imagePtr + i5);
          dotProduct += *(kernelPtr + k6) * *(imagePtr + i6);
          dotProduct += *(kernelPtr + k7) * *(imagePtr + i7);
          dotProduct += *(kernelPtr + k8) * *(imagePtr + i8);
          resultImage(resultIndex) = dotProduct;
          ++imagePtr;
          ++resultIndex;
        }
      }
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::Stencil2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Stencil2DTest currentTest;

}

#endif
