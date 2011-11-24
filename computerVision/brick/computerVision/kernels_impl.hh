/**
***************************************************************************
* @file brick/computerVision/kernels.hh
*
* Header file declaring some stock kernels to be used in client code.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KERNELS_IMPL_HH
#define BRICK_COMPUTERVISION_KERNELS_IMPL_HH

// This file is included by kernels.hh, and should not be directly included
// by user code, so no need to include kernel.hh here.
// 
// #include <brick/computerVision/kernels.hh>

#include <cmath>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace computerVision {

    /// @cond privateCode
    namespace privateCode {
      
      template <class TYPE>
      brick::numeric::Array1D<TYPE>
      getGaussian1D(size_t size, double sigma, bool normalize=false)
      {
        const double myPi = 3.14159265359;
        brick::numeric::Array1D<TYPE> result(size);
        double x = (1.0 - static_cast<double>(size))/2.0;
        double twoSigmaSq = 2.0 * sigma * sigma;
        double k = 1.0 / (std::sqrt(2.0 * myPi) * sigma);
        for(size_t index0 = 0; index0 < size; ++index0) {
          result[index0] = k * exp(-x * x / twoSigmaSq);
          x += 1.0;
        }
        if(normalize) {
          result /= brick::numeric::sum<TYPE>(result);
        }
        return result;
      }
      
    } // namespace privateCode
    /// @endcond


    // This function generates and returns a separable Gaussian kernel.
    template<class KERNEL_TYPE>
    Kernel<KERNEL_TYPE>
    getGaussianKernel(double rowSigma, double columnSigma)
    {
      size_t rows = static_cast<size_t>(6.0 * rowSigma + 1.0);
      size_t columns = static_cast<size_t>(6.0 * columnSigma + 1.0);
      if(rows % 2 == 0) {
        ++rows;
      }
      if(columns % 2 == 0) {
        ++columns;
      }
      return getGaussianKernel<KERNEL_TYPE>(
        rows, columns, rowSigma, columnSigma);
    }

    
    // This function generates and returns a separable Gaussian kernel.
    template<class KERNEL_TYPE>
    Kernel<KERNEL_TYPE>
    getGaussianKernel(size_t rows, size_t columns,
                      double rowSigma, double columnSigma)
    {
      // Argument checking.
      if(rows == 0 || columns == 0) {
        BRICK_THROW(brick::common::ValueException, "getGaussianKernel()",
                  "Arguments rows and columns may not have value of zero.");
      }
      if(rowSigma < 0.0) {
        rowSigma = rows / 6.0;
      }
      if(columnSigma < 0.0) {
        columnSigma = columns / 6.0;
      }

      brick::numeric::Array1D<KERNEL_TYPE> rowComponent =
        privateCode::getGaussian1D<KERNEL_TYPE>(columns, columnSigma, true);
      brick::numeric::Array1D<KERNEL_TYPE> columnComponent =
        privateCode::getGaussian1D<KERNEL_TYPE>(rows, rowSigma, true);
      return Kernel<KERNEL_TYPE>(rowComponent, columnComponent);
    }


  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KERNELS_IMPL_HH */
