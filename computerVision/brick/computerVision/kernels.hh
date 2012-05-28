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

#ifndef BRICK_COMPUTERVISION_KERNELS_HH
#define BRICK_COMPUTERVISION_KERNELS_HH

#include <brick/computerVision/kernel.hh>

namespace brick {

  namespace computerVision {
    
    /** 
     * This function generates and returns a separable Gaussian kernel
     * of user-specified standard deviation.  The kernel will be sized
     * so that the number of rows is equal to the smallest odd number
     * greater than 6.0 * rowSigma, and the number of columns is equal
     * to the smallest odd number greater than 6.0 * columnSigma.  The
     * kernel will be normalized prior to return, so that the sum of
     * its elements is equal to one.
     *
     * Use this function like this:
     *
     * @code
     *   Kernel<double> kernel = getGaussianKernel<double>(1.0, 1.0);
     * @endcode
     * 
     * @param rowSigma This argument specifies the standard deviation
     * of the kernel in the Y direction.
     *
     * @param columnSigma This argument specifies the standard deviation of
     * the kernel in the X direction.
     * 
     * @param normalize This argument indicates whether the resulting
     * gaussian should be normalized so that its integral has known
     * value.  Note that, as a seperable kernel, row and column
     * components will be normalized separately, so that the integral
     * of the effective 2D kernel will be (rowNormalizationTarget *
     * columnNormalizationTarget).
     * 
     * @param rowNormalizationTarget If the row component of the
     * kernel is normalized, it will be normalize so that its integral
     * is as close as possible to this value.
     * 
     * @param columnNormalizationTarget If the column component of the
     * kernel is normalized, it will be normalize so that its integral
     * is as close as possible to this value.
     * 
     * @return The return value is the kernel itself.
     */
    template<class KERNEL_TYPE>
    Kernel<KERNEL_TYPE>
    getGaussianKernel(double rowSigma, double columnSigma,
                      bool normalize = true,
                      KERNEL_TYPE rowNormalizationTarget = 1,
                      KERNEL_TYPE columnNormalizationTarget = 1);
                      


    /** 
     * This function generates and returns a separable Gaussian kernel
     * with the same variance along each axis.  The kernel will be
     * normalized prior to return, so that the sum of its elements is
     * equal to one.
     *
     * Use this function like this:
     *
     * @code
     *   Kernel<double> kernel = getGaussianKernel<double>(5, 5);
     * @endcode
     * 
     * @param rows This argument specifies how many rows the kernel
     * should have.
     * 
     * @param columns This argument specifies how many columns the
     * kernel should have.
     * 
     * @param rowSigma This argument specifies the standard deviation
     * of the kernel in the Y direction.  If sigma is less than 0.0,
     * then it will be automatically reset to rows / 6.0.
     *
     * @param columnSigma This argument specifies the standard deviation of
     * the kernel in the X direction.  If sigma is less than 0.0, then
     * it will be automatically reset columns / 6.0.
     * 
     * @param normalize This argument indicates whether the resulting
     * gaussian should be normalized so that its integral has known
     * value.  Note that, as a seperable kernel, row and column
     * components will be normalized separately, so that the integral
     * of the effective 2D kernel will be (rowNormalizationTarget *
     * columnNormalizationTarget).
     * 
     * @param rowNormalizationTarget If the row component of the
     * kernel is normalized, it will be normalize so that its integral
     * is as close as possible to this value.
     * 
     * @param columnNormalizationTarget If the column component of the
     * kernel is normalized, it will be normalize so that its integral
     * is as close as possible to this value.
     * 
     * @return The return value is the kernel itself.
     */
    template<class KERNEL_TYPE>
    Kernel<KERNEL_TYPE>
    getGaussianKernelBySize(size_t rows, size_t columns,
                      double rowSigma=-1.0, double columnSigma=-1.0,
                      bool normalize = true,
                      KERNEL_TYPE rowNormalizationTarget = 1,
                      KERNEL_TYPE columnNormalizationTarget = 1);

  } // namespace computerVision

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/kernels_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KERNELS_HH */
