/**
***************************************************************************
* @file brick/computerVision/keypointSelectorHarris_impl.hh
*
* Header file defining a class template for selecting stable keypoints
* from an image.
*
* Copyright (C) 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORHARRIS_IMPL_HH
#define BRICK_COMPUTERVISION_KEYPOINTSELECTORHARRIS_IMPL_HH

// This file is included by keypointSelectorHarris.hh, and should not be directly included
// by user code, so no need to include keypointSelectorHarris.hh here.
//
// #include <brick/computerVision/keypointSelectorHarris.hh>

#include <brick/common/mathFunctions.hh>
#include <brick/computerVision/imageFilter.hh>
#include <brick/numeric/bilinearInterpolator.hh>
#include <brick/numeric/filter.hh>
#include <brick/numeric/subpixelInterpolate.hh>

#ifndef BRICK_COMPUTERVISION_HARRIS_PEDANTIC
#define BRICK_COMPUTERVISION_HARRIS_PEDANTIC 0
#endif /* #ifndef BRICK_COMPUTERVISION_HARRIS_PEDANTIC */

namespace brick {

  namespace computerVision {


    template <class CoordinateType>
    template <class FloatType>
    void
    KeypointHarris<CoordinateType>::
    getCovariance(FloatType& c00, FloatType& c01, FloatType& c11)
    {
      brick::common::Int32 determinant = m_xx * m_yy - m_xy * m_xy;
      c00 = FloatType(m_yy) / determinant;
      c01 = -FloatType(m_xy) / determinant;
      c11 = FloatType(m_yy) / determinant;
    }


    template <class FloatType>
    KeypointSelectorHarris<FloatType>::
    KeypointSelectorHarris(FloatType kappa,
                           FloatType sigma)
      : m_harrisIndicators(),
        m_gradientXX(),
        m_gradientXY(),
        m_gradientYY(),
        m_searchRegionCorner0(0, 0),
        m_searchRegionCorner1(0, 0),
        m_kappa(kappa),
        m_sigma(sigma)
    {
      // Empty.
    }


    // Return the keypoints detected during the most recent call to
    // member function setImage().
    template <class FloatType>
    std::vector< KeypointHarris<brick::common::Int32> >
    KeypointSelectorHarris<FloatType>::
    getKeypoints() const
    {
      std::vector< KeypointHarris<brick::common::Int32> > keypointVector;

      // Keep only those pixels that are bigger than their neighbors.
      this->getKeypoints(std::back_inserter(keypointVector));
      return keypointVector;
    }


    // Return the keypoints detected during the most recent call to
    // member function setImage().
    template <class FloatType>
    template <class Iter>
    void
    KeypointSelectorHarris<FloatType>::
    getKeypoints(Iter iterator, FloatType /* threshold */) const
    {
      // Get bounds on the region of valid pixels for this calculation.
      unsigned int const startRow    = m_searchRegionCorner0.getRow();
      unsigned int const stopRow     = m_searchRegionCorner1.getRow();
      unsigned int const startColumn = m_searchRegionCorner0.getColumn();
      unsigned int const stopColumn  = m_searchRegionCorner1.getColumn();

      // Iterate over all pixels in the valid region.
      unsigned int rowStep = m_harrisIndicators.getRowStep();
      for(unsigned int row = stopRow - 1; row >= startRow; --row) {
        brick::numeric::Array1D<AccumulatedType> xxRow =
          m_gradientXX.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xyRow =
          m_gradientXY.getRow(row);
        brick::numeric::Array1D<AccumulatedType> yyRow =
          m_gradientYY.getRow(row);
        brick::numeric::Array1D<FloatType> harrisRow =
          m_harrisIndicators.getRow(row);
        for(unsigned int column = stopColumn - 1; column >= startColumn;
            --column) {
          FloatType* candidatePtr = &(harrisRow[column]);
          if(// (*candidatePtr > threshold) &&
            (*candidatePtr > *(candidatePtr + 1))
             && (*candidatePtr > *(candidatePtr - 1))
             && (*candidatePtr > *(candidatePtr + rowStep))
             && (*candidatePtr > *(candidatePtr - rowStep))
             && (*candidatePtr > *(candidatePtr + rowStep + 1))
             && (*candidatePtr > *(candidatePtr - rowStep - 1))
             && (*candidatePtr > *(candidatePtr + rowStep - 1))
             && (*candidatePtr > *(candidatePtr - rowStep + 1))) {
            *(iterator++) = KeypointHarris<brick::common::Int32>(
              row, column, *candidatePtr,
              xxRow[column], yyRow[column], xyRow[column]);
          }
        }
      }
    }


    template <class FloatType>
    std::vector< KeypointHarris<FloatType> >
    KeypointSelectorHarris<FloatType>::
    getKeypointsGeneralPosition() const
    {
      std::vector< KeypointHarris<FloatType> > keypointVector;

      // Keep only those pixels that are bigger than their neighbors.
      this->getKeypointsGeneralPosition(std::back_inserter(keypointVector));
      return keypointVector;
    }


    // Return the keypoints detected during the most recent call to
    // member function setImage().
    template <class FloatType>
    template <class Iter>
    void
    KeypointSelectorHarris<FloatType>::
    getKeypointsGeneralPosition(Iter iterator) const
    {
      // Sanity check region to make sure shrinking the search region
      // isn't going to break anything.
      if(0 == m_searchRegionCorner1.getRow()
         || 0 == m_searchRegionCorner1.getColumn()) {
        return;
      }

      // Get bounds on the region of valid pixels for this
      // calculation.  Shrink search region by one pixel to allow
      // non-max suppression.
      unsigned int const startRow    = m_searchRegionCorner0.getRow() + 1;
      unsigned int const stopRow     = m_searchRegionCorner1.getRow() - 1;
      unsigned int const startColumn = m_searchRegionCorner0.getColumn() + 1;
      unsigned int const stopColumn  = m_searchRegionCorner1.getColumn() - 1;

      // These interpolators will make things easy on us later when
      // trying to gather information about keypoints in non-integral
      // positions.
      brick::numeric::BilinearInterpolator<AccumulatedType, FloatType>
        xxInterpolator(m_gradientXX);
      brick::numeric::BilinearInterpolator<AccumulatedType, FloatType>
        xyInterpolator(m_gradientXY);
      brick::numeric::BilinearInterpolator<AccumulatedType, FloatType>
        yyInterpolator(m_gradientYY);

      // Iterate over all pixels in the valid region.
      unsigned int rowStep = m_harrisIndicators.getRowStep();
      for(unsigned int row = stopRow - 1; row >= startRow; --row) {
        brick::numeric::Array1D<AccumulatedType> xxRow =
          m_gradientXX.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xyRow =
          m_gradientXY.getRow(row);
        brick::numeric::Array1D<AccumulatedType> yyRow =
          m_gradientYY.getRow(row);
        brick::numeric::Array1D<FloatType> harrisRow =
          m_harrisIndicators.getRow(row);
        for(unsigned int column = stopColumn - 1; column >= startColumn;
            --column) {
          FloatType* candidatePtr = &(harrisRow[column]);
          if((*candidatePtr > *(candidatePtr + 1))
             && (*candidatePtr > *(candidatePtr - 1))
             && (*candidatePtr > *(candidatePtr + rowStep))
             && (*candidatePtr > *(candidatePtr - rowStep))
             && (*candidatePtr > *(candidatePtr + rowStep + 1))
             && (*candidatePtr > *(candidatePtr - rowStep - 1))
             && (*candidatePtr > *(candidatePtr + rowStep - 1))
             && (*candidatePtr > *(candidatePtr - rowStep + 1))) {

            FloatType rowCoordinate;
            FloatType columnCoordinate;
            FloatType extremeValue;
            if(brick::numeric::subpixelInterpolate(
                 FloatType(row), FloatType(column),
                 *(candidatePtr - rowStep - 1), *(candidatePtr - rowStep),
                 *(candidatePtr - rowStep + 1),
                 *(candidatePtr - 1), *candidatePtr, *(candidatePtr + 1),
                 *(candidatePtr + rowStep - 1), *(candidatePtr + rowStep),
                 *(candidatePtr + rowStep + 1),
                 rowCoordinate, columnCoordinate, extremeValue)) {

              // Sanity check interpolation result, as subpixelInterpolate()
              // does not do so.
              if(common::absoluteValue(rowCoordinate - row) < 1.0
                 && common::absoluteValue(columnCoordinate - column) < 1.0) {

                *(iterator++) = KeypointHarris<FloatType>(
                  rowCoordinate, columnCoordinate, extremeValue,
                  xxInterpolator(rowCoordinate, columnCoordinate),
                  xyInterpolator(rowCoordinate, columnCoordinate),
                  yyInterpolator(rowCoordinate, columnCoordinate));

              }
            }
          }
        }
      }
    }


    template <class FloatType>
    void
    KeypointSelectorHarris<FloatType>::
    setImage(Image<GRAY8> const& inImage)
    {
      // Start by blurring the image slightly.  This should be
      // optional in future implementations.  This call makes an
      // integer valued approximation to a 2D Gaussian kernel,
      // normalized so that it integrates to approximately 256 * 256 =
      // 65536.  The size of this filter affects the region of the
      // image for which the Harris corner metric will be valid, so we
      // update the search region here as well.
      Kernel<brick::common::Int32> gaussian =
        getGaussianKernelBySize<brick::common::Int32>(
          size_t(5), size_t(5), -1.0, -1.0, true, 256, 256);
      m_searchRegionCorner0.setValue(2, 2);
      m_searchRegionCorner1.setValue(inImage.rows() - 2, inImage.columns() - 2);

      // Here we do the convolution.  Fortunately, 2D Gaussians are
      // separable, so we can do the row and column convolutions
      // separately.  In a more cache-coherent implementation, we'd
      // combine these operations.
      brick::numeric::Array2D<brick::common::Int32> workspace(
        inImage.rows(), inImage.columns());
      workspace = 0;
      Image<GRAY_SIGNED32> blurredImage(inImage.rows(), inImage.columns());
      brick::numeric::filterRows(
        workspace, inImage, gaussian.getRowComponent());
      brick::numeric::filterColumns(
        blurredImage, workspace, gaussian.getColumnComponent());

      // The convolution has increased the average value of the image
      // by a factor of about 65536.  If we proceed without scaling
      // down, we'll overflow our 32 bit ints!  Fix this here.
      blurredImage >>= 16;

      // Now compute products of gradients at each pixel.  Note that
      // the gradient calculation is not valid at the very edge of the
      // window (where the gradient filter extends into no-mans-land)
      // so search region will be reduced slightly by this call.
      this->computeGradients(
        blurredImage, m_gradientXX, m_gradientXY, m_gradientYY,
        m_searchRegionCorner0, m_searchRegionCorner1);

      // To get rotation independent Harris corners, the integration
      // for each pixel needs to be weighted with a circularly
      // symmetric window.  We'll use another gaussian blur.  Remember
      // that we're convolving with the squares of gradients.  To
      // avoid risk of overflow, this kernel is normalized so it
      // integrates to a smaller value (45 * 45 = 2025).
      gaussian = getGaussianKernelBySize<brick::common::Int32>(
        size_t(11), size_t(11), -1.0, -1.0, true, 45, 45);
      unsigned int radius = 5;

      // Do the weighted integration by means of separable convolution.
      brick::numeric::Array2D<AccumulatedType> proxyArray0;
      brick::numeric::Array2D<AccumulatedType> proxyArray1;

      // We need these proxies because otherwise g++ complains about
      // passing an rvalue to initialize non-const reference in the
      // filter*() calls.
      proxyArray0 = workspace.getRegion(m_searchRegionCorner0,
                                        m_searchRegionCorner1);
      proxyArray1 = m_gradientXX.getRegion(m_searchRegionCorner0,
                                           m_searchRegionCorner1),
      brick::numeric::filterRows(
        proxyArray0, proxyArray1, gaussian.getRowComponent());
      brick::numeric::filterColumns(
        proxyArray1, proxyArray0, gaussian.getColumnComponent());

      proxyArray1 = m_gradientXY.getRegion(m_searchRegionCorner0,
                                           m_searchRegionCorner1),
      brick::numeric::filterRows(
        proxyArray0, proxyArray1, gaussian.getRowComponent());
      brick::numeric::filterColumns(
        proxyArray1, proxyArray0, gaussian.getColumnComponent());

      proxyArray1 = m_gradientYY.getRegion(m_searchRegionCorner0,
                                           m_searchRegionCorner1),
      brick::numeric::filterRows(
        proxyArray0, proxyArray1, gaussian.getRowComponent());
      brick::numeric::filterColumns(
        proxyArray1, proxyArray0, gaussian.getColumnComponent());

      // The size of this filter affects the region of the image for
      // which the Harris corner metric will be valid, so we update
      // the search region here as well.
      m_searchRegionCorner0.setValue(
        m_searchRegionCorner0.getRow() + radius,
        m_searchRegionCorner0.getColumn() + radius);
      m_searchRegionCorner1.setValue(
        m_searchRegionCorner1.getRow() - radius,
        m_searchRegionCorner1.getColumn() - radius);

      // Again we have to worry about overflow.  Rather than guessing,
      // just search to find the max value, and rescale based on that.
      // Ideally, this search would be done during the convolution
      // above.
      brick::common::Int32 maxVal =
        brick::numeric::maximum(
          m_gradientXX.getRegion(m_searchRegionCorner0,
                                 m_searchRegionCorner1));
      maxVal = std::max(
        maxVal, brick::numeric::maximum(
          m_gradientXY.getRegion(m_searchRegionCorner0,
                                 m_searchRegionCorner1)));
      maxVal = std::max(
        maxVal, brick::numeric::maximum(
          m_gradientYY.getRegion(m_searchRegionCorner0,
                                 m_searchRegionCorner1)));
      brick::common::Int32 divisor = maxVal / 65535 + 1;

      // Rescale each of the images.  The exact rescaling factor
      // doesn't matter too much here, as all of the processing
      // downstream from here is theoretically scale-independent.  The
      // bigger divisor is, however, the more precision we lose, so we
      // want divisor to be small.
      m_gradientXX.getRegion(m_searchRegionCorner0, m_searchRegionCorner1)
        /= divisor;
      m_gradientXY.getRegion(m_searchRegionCorner0, m_searchRegionCorner1)
        /= divisor;
      m_gradientYY.getRegion(m_searchRegionCorner0, m_searchRegionCorner1)
        /= divisor;

      // Compute the Harris indicator for each pixel in the region.
      this->computeHarrisIndicators(
        m_gradientXX, m_gradientXY, m_gradientYY, m_harrisIndicators,
        m_searchRegionCorner0, m_searchRegionCorner1);
    }


    // ============== Private member functions below this line ==============

    template <class FloatType>
    void
    KeypointSelectorHarris<FloatType>::
    computeGradients(
      Image<GRAY_SIGNED32> const& inImage,
      brick::numeric::Array2D<AccumulatedType>& gradientXX,
      brick::numeric::Array2D<AccumulatedType>& gradientXY,
      brick::numeric::Array2D<AccumulatedType>& gradientYY,
      brick::numeric::Index2D& corner0,
      brick::numeric::Index2D& corner1)
    {
      int const rowStep     = inImage.getRowStep();
      int const startRow    = corner0.getRow();
      int const stopRow     = corner1.getRow();
      int const startColumn = corner0.getColumn();
      int const stopColumn  = corner1.getColumn();

      // Adjust array sizes, if necessary.
      if((gradientXX.rows() != inImage.rows())
         || (gradientXX.columns() != inImage.columns())) {
        gradientXX.reinit(inImage.rows(), inImage.columns());
        gradientXY.reinit(inImage.rows(), inImage.columns());
        gradientYY.reinit(inImage.rows(), inImage.columns());
      }

#if BRICK_COMPUTERVISION_HARRIS_PEDANTIC
      // Zero out untouched (and unused) pixels.
      gradientXX = 0;
      gradientXY = 0;
      gradientYY = 0;
#endif /* #if BRICK_COMPUTERVISION_HARRIS_PEDANTIC */

      for(int row = stopRow - 1; row >= startRow; --row) {
        brick::numeric::Array1D<brick::common::Int32> imageRow =
          inImage.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xxRow =
          gradientXX.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xyRow =
          gradientXY.getRow(row);
        brick::numeric::Array1D<AccumulatedType> yyRow =
          gradientYY.getRow(row);

        // Note that in the filter code below, we actually index
        // _outside_ of the data pointed to by imageRow.  This is "ok"
        // because we know that imageRow points into inImage (and we
        // know the memory layout of inImage), so we're sure there's
        // valid data, even past the end or beginning of the row.
        for(int column = stopColumn - 1; column >= startColumn;
            --column) {
          brick::common::Int32 gradientX =
            (((static_cast<brick::common::Int32>(imageRow[column + 1])
               - static_cast<brick::common::Int32>(imageRow[column - 1])) << 1)
             + (static_cast<brick::common::Int32>(
                  *(imageRow.getData(column) + rowStep + 1))
                - static_cast<brick::common::Int32>(
                  *(imageRow.getData(column) + rowStep - 1)))
             + (static_cast<brick::common::Int32>(
                  *(imageRow.getData(column) - rowStep + 1))
                - static_cast<brick::common::Int32>(
                  *(imageRow.getData(column) - rowStep - 1)))); /* >> 2; */
          brick::common::Int32 gradientY =
            (((static_cast<brick::common::Int32>(
                 *(imageRow.getData(column) + rowStep))
               - static_cast<brick::common::Int32>(
                 *(imageRow.getData(column) - rowStep))) << 1)
             + (static_cast<brick::common::Int32>(
                  *(imageRow.getData(column) + rowStep - 1))
                - static_cast<brick::common::Int32>(
                  *(imageRow.getData(column) - rowStep - 1)))
             + (static_cast<brick::common::Int32>(
                  *(imageRow.getData(column) + rowStep + 1))
                - static_cast<brick::common::Int32>(
                  *(imageRow.getData(column) - rowStep + 1)))); /* >> 2 */;
          xxRow[column] = (gradientX * gradientX); /* >> 4; */
          xyRow[column] = (gradientX * gradientY); /* >> 4; */
          yyRow[column] = (gradientY * gradientY); /* >> 4; */
        }
      }

      // Adjust search region to reflect the fact that the gradient
      // result is not valid on border pixels.
      corner0.setValue(startRow + 1, startColumn + 1);
      corner1.setValue(stopRow  - 1, stopColumn  - 1);
    }


    template <class FloatType>
    void
    KeypointSelectorHarris<FloatType>::
    computeHarrisIndicators(
      brick::numeric::Array2D<AccumulatedType> const& gradientXX,
      brick::numeric::Array2D<AccumulatedType> const& gradientXY,
      brick::numeric::Array2D<AccumulatedType> const& gradientYY,
      brick::numeric::Array2D<FloatType>& harrisIndicators,
      brick::numeric::Index2D const& corner0,
      brick::numeric::Index2D const& corner1)
    {
      // Get bounds on the region of valid pixels for this calculation.
      unsigned int const startRow    = corner0.getRow();
      unsigned int const stopRow     = corner1.getRow();
      unsigned int const startColumn = corner0.getColumn();
      unsigned int const stopColumn  = corner1.getColumn();

      // Adjust array size, if necessary.
      if((harrisIndicators.rows() != gradientXX.rows())
         || (harrisIndicators.columns() != gradientXX.columns())) {
        harrisIndicators.reinit(gradientXX.rows(), gradientXX.columns());
      }

      // Iterate over all pixels in the valid region.
      for(unsigned int row = stopRow - 1; row >= startRow; --row) {
        brick::numeric::Array1D<AccumulatedType> xxRow =
          gradientXX.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xyRow =
          gradientXY.getRow(row);
        brick::numeric::Array1D<AccumulatedType> yyRow =
          gradientYY.getRow(row);
        brick::numeric::Array1D<FloatType> harrisRow =
          harrisIndicators.getRow(row);
        for(unsigned int column = stopColumn - 1; column >= startColumn;
            --column) {
          AccumulatedType determinant = (xxRow[column] * yyRow[column]
                                - xyRow[column] * xyRow[column]);
          AccumulatedType trace = xxRow[column] + yyRow[column];
          harrisRow[column] = determinant - m_kappa * trace * trace;
        }
      }
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORHARRIS_IMPL_HH */
