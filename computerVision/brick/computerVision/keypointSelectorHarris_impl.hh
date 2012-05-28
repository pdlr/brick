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

#include <brick/computerVision/imageFilter.hh>
#include <brick/numeric/bilinearInterpolator.hh>
#include <brick/numeric/subpixelInterpolate.hh>

#ifndef BRICK_COMPUTERVISION_HARRIS_PEDANTIC
#define BRICK_COMPUTERVISION_HARRIS_PEDANTIC 1
#endif /* #ifndef BRICK_COMPUTERVISION_HARRIS_PEDANTIC */

#define BRICK_COMPUTERVISION_HARRIS_DEBUG 1

#if BRICK_COMPUTERVISION_HARRIS_DEBUG
#include <brick/computerVision/sobel.hh>
#include <brick/computerVision/utilities.hh>
#include <brick/utilities/timeUtilities.hh>
#endif


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
      unsigned int rowStep = m_harrisIndicators.getRowStep();
      for(unsigned int row = m_harrisIndicators.rows() - 2; row > 1; --row) {
        brick::numeric::Array1D<AccumulatedType> xxRow =
          m_gradientXX.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xyRow =
          m_gradientXY.getRow(row);
        brick::numeric::Array1D<AccumulatedType> yyRow =
          m_gradientYY.getRow(row);
        brick::numeric::Array1D<FloatType> harrisRow =
          m_harrisIndicators.getRow(row);
        for(unsigned int column = m_harrisIndicators.columns() - 2; column > 1;
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
      brick::numeric::BilinearInterpolator<AccumulatedType, FloatType>
        xxInterpolator(m_gradientXX);
      brick::numeric::BilinearInterpolator<AccumulatedType, FloatType>
        xyInterpolator(m_gradientXY);
      brick::numeric::BilinearInterpolator<AccumulatedType, FloatType>
        yyInterpolator(m_gradientYY);
      
      unsigned int rowStep = m_harrisIndicators.getRowStep();
      for(unsigned int row = m_harrisIndicators.rows() - 2; row > 1; --row) {
        brick::numeric::Array1D<AccumulatedType> xxRow =
          m_gradientXX.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xyRow =
          m_gradientXY.getRow(row);
        brick::numeric::Array1D<AccumulatedType> yyRow =
          m_gradientYY.getRow(row);
        brick::numeric::Array1D<FloatType> harrisRow =
          m_harrisIndicators.getRow(row);
        for(unsigned int column = m_harrisIndicators.columns() - 2; column > 1;
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


    template <class FloatType>
    void
    KeypointSelectorHarris<FloatType>::
    setImage(Image<GRAY8> const& inImage)
    {
      Kernel<brick::common::Int32> gaussian = 
        getGaussianKernelBySize<brick::common::Int32>(
          size_t(5), size_t(5), -1.0, -1.0, true, 256, 256);
      Image<GRAY_SIGNED32> blurredImage = filter2D<GRAY_SIGNED32>(
        gaussian, inImage, 0.0);
      blurredImage /= 65536;
      this->computeGradients(
        blurredImage, m_gradientXX, m_gradientXY, m_gradientYY);


#if 1 
      Kernel<brick::common::Int32> gaussian2_8 =
        getGaussianKernelBySize<brick::common::Int32>(
          size_t(11), size_t(11), -1.0, -1.0, true, 45, 45);

      m_gradientXX = filter2D<GRAY_SIGNED32, GRAY_SIGNED32>(
        gaussian2_8, Image<GRAY_SIGNED32>(m_gradientXX), 0);
      m_gradientXY = filter2D<GRAY_SIGNED32, GRAY_SIGNED32>(
        gaussian2_8, Image<GRAY_SIGNED32>(m_gradientXY), 0);
      m_gradientYY = filter2D<GRAY_SIGNED32, GRAY_SIGNED32>(
        gaussian2_8, Image<GRAY_SIGNED32>(m_gradientYY), 0);

      brick::common::Int32 maxVal =
        brick::numeric::maximum(m_gradientXX.ravel());
      maxVal = std::max(
        maxVal, brick::numeric::maximum(m_gradientXY.ravel()));
      maxVal = std::max(
        maxVal, brick::numeric::maximum(m_gradientYY.ravel()));
      std::cout << "========" << maxVal << "========" << std::endl;
      brick::common::Int32 divisor = maxVal / 65535 + 1;
      m_gradientXX /= divisor;
      m_gradientXY /= divisor;
      m_gradientYY /= divisor;
#else      
      double t0 = utilities::getCurrentTime();

      // Apply a window function to the products of gradients.
      Image<GRAY_SIGNED32> xxImage(m_gradientXX);
      Image<GRAY_SIGNED32> xyImage(m_gradientXY);
      Image<GRAY_SIGNED32> yyImage(m_gradientYY);
      Image<GRAY_SIGNED32> workingImage(inImage.rows(), inImage.columns());
      filterRowsBinomial<brick::common::Int32>(
        workingImage, xxImage,      this->m_sigma);
      filterColumnsBinomial<brick::common::Int32>(
        xxImage,      workingImage, this->m_sigma);
      filterRowsBinomial<brick::common::Int32>(
        workingImage, xyImage,      this->m_sigma);
      filterColumnsBinomial<brick::common::Int32>(
        xyImage,      workingImage, this->m_sigma);
      filterRowsBinomial<brick::common::Int32>(
        workingImage, yyImage,      this->m_sigma);
      filterColumnsBinomial<brick::common::Int32>(
        yyImage,      workingImage, this->m_sigma);
      double t1 = utilities::getCurrentTime();
      std::cout << "ET: " << t1 - t0 << std::endl;
#endif

      // Compute the Harris indicator for each pixel in the region.
      this->computeHarrisIndicators(m_gradientXX, m_gradientXY, m_gradientYY,
                                    m_harrisIndicators);
    }


    

    // ============== Private member functions below this line ==============

    template <class FloatType>
    void
    KeypointSelectorHarris<FloatType>::
    computeGradients(
      Image<GRAY_SIGNED32> const& inImage,
      brick::numeric::Array2D<AccumulatedType>& gradientXX,
      brick::numeric::Array2D<AccumulatedType>& gradientXY,
      brick::numeric::Array2D<AccumulatedType>& gradientYY)
    {
      unsigned int const rowStep = inImage.getRowStep();

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
      
      for(unsigned int row = inImage.rows() - 2; row > 1; --row) {
        brick::numeric::Array1D<brick::common::Int32> imageRow =
          inImage.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xxRow =
          gradientXX.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xyRow =
          gradientXY.getRow(row);
        brick::numeric::Array1D<AccumulatedType> yyRow =
          gradientYY.getRow(row);
        for(unsigned int column = inImage.columns() - 2; column > 1; --column) {
          brick::common::Int32 gradientX =
            (((static_cast<brick::common::Int32>(imageRow[column + 1])
               - static_cast<brick::common::Int32>(imageRow[column - 1])) << 1)
             + (static_cast<brick::common::Int32>(
                  imageRow[column + rowStep + 1])
                - static_cast<brick::common::Int32>(
                  imageRow[column + rowStep - 1]))
            + (static_cast<brick::common::Int32>(
                 imageRow[column - rowStep + 1])
               - static_cast<brick::common::Int32>(
                 imageRow[column - rowStep - 1]))); /* >> 2; */
          brick::common::Int32 gradientY =
            (((static_cast<brick::common::Int32>(
                 imageRow[column + rowStep])
               - static_cast<brick::common::Int32>(
                 imageRow[column - rowStep])) << 1)
             + (static_cast<brick::common::Int32>(
                  imageRow[column + rowStep - 1])
                - static_cast<brick::common::Int32>(
                  imageRow[column - rowStep - 1]))
             + (static_cast<brick::common::Int32>(
                  imageRow[column + rowStep + 1])
                - static_cast<brick::common::Int32>(
                  imageRow[column - rowStep + 1]))); /* >> 2 */;
          xxRow[column] = (gradientX * gradientX); /* >> 4; */
          xyRow[column] = (gradientX * gradientY); /* >> 4; */
          yyRow[column] = (gradientY * gradientY); /* >> 4; */
        }
      }
    }
    
    
    template <class FloatType>
    void
    KeypointSelectorHarris<FloatType>::
    computeHarrisIndicators(
      brick::numeric::Array2D<AccumulatedType> const& gradientXX,
      brick::numeric::Array2D<AccumulatedType> const& gradientXY,
      brick::numeric::Array2D<AccumulatedType> const& gradientYY,
      brick::numeric::Array2D<FloatType>& harrisIndicators)
    {
      // Adjust array size, if necessary.
      if((harrisIndicators.rows() != gradientXX.rows())
         || (harrisIndicators.columns() != gradientXX.columns())) {
        harrisIndicators.reinit(gradientXX.rows(), gradientXX.columns());
      }

      for(unsigned int row = harrisIndicators.rows() - 2; row > 1; --row) {
        brick::numeric::Array1D<AccumulatedType> xxRow =
          gradientXX.getRow(row);
        brick::numeric::Array1D<AccumulatedType> xyRow =
          gradientXY.getRow(row);
        brick::numeric::Array1D<AccumulatedType> yyRow =
          gradientYY.getRow(row);
        brick::numeric::Array1D<FloatType> harrisRow =
          harrisIndicators.getRow(row);
        for(unsigned int column = harrisIndicators.columns() - 2; column > 1;
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
