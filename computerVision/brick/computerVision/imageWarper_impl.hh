/**
***************************************************************************
* @file brick/computerVision/imageWarper.hh
*
* Header file defining inline and template functions declared in
* imageWarper.hh.
*
* Copyright (C) 2009,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEWARPER_IMPL_HH
#define BRICK_COMPUTERVISION_IMAGEWARPER_IMPL_HH

// This file is included by imageWarper.hh, and should not be directly included
// by user code, so no need to include imageWarper.hh here.
//
// #include <brick/computerVision/imageWarper.hh>

#include <brick/numeric/vector2D.hh>
#include <brick/common/mathFunctions.hh>

namespace brick {

  namespace computerVision {


    template<class NumericType, class TransformFunctor>
    ImageWarper<NumericType, TransformFunctor>::
    ImageWarper()
      : m_inputColumns(0),
        m_inputRows(0),
        m_lookupTable()
    {
      // Empty.
    }


    template<class NumericType, class TransformFunctor>
    ImageWarper<NumericType, TransformFunctor>::
    ImageWarper(size_t inputRows, size_t inputColumns,
                size_t outputRows, size_t outputColumns,
                TransformFunctor transformer)
      : m_inputColumns(inputColumns),
        m_inputRows(inputRows),
        m_lookupTable(outputRows, outputColumns)
    {
      for(size_t row = 0; row < outputRows; ++row) {
        for(size_t column = 0; column < outputColumns; ++column) {
          brick::numeric::Vector2D<NumericType> outputCoord(column, row);
          brick::numeric::Vector2D<NumericType> inputCoord =
            transformer(outputCoord);
          SampleInfo& sampleInfo = m_lookupTable(row, column);
          if((inputCoord.y() >= 0.0)
             && (inputCoord.x() >= 0.0)
             && (inputCoord.y() < static_cast<NumericType>(inputRows - 1))
             && (inputCoord.x() < static_cast<NumericType>(inputColumns - 1))) {

            NumericType intPart;
            NumericType xFrac;
            NumericType yFrac;

            brick::common::splitFraction(inputCoord.x(), intPart, xFrac);
            size_t i0 = static_cast<size_t>(intPart);

            brick::common::splitFraction(inputCoord.y(), intPart, yFrac);
            size_t j0 = static_cast<size_t>(intPart);

            NumericType oneMinusXFrac = 1.0 - xFrac;
            NumericType oneMinusYFrac = 1.0 - yFrac;
            sampleInfo.c00 = oneMinusXFrac * oneMinusYFrac;
            sampleInfo.c01 = xFrac * oneMinusYFrac;
            sampleInfo.c10 = oneMinusXFrac * yFrac;
            sampleInfo.c11 = xFrac * yFrac;
            sampleInfo.index00 = inputColumns * j0 + i0;
            sampleInfo.isInBounds = true;
          } else {
            sampleInfo.isInBounds = false;
          }
        }
      }
    }


    // Destroys the ImageWarper instance and deletes the internal data
    // store.
    template<class NumericType, class TransformFunctor>
    ImageWarper<NumericType, TransformFunctor>::
    ~ImageWarper()
    {
      // Empty.
    }



    template<class NumericType, class TransformFunctor>
    template <ImageFormat InputFormat, ImageFormat OutputFormat>
    Image<OutputFormat>
    ImageWarper<NumericType, TransformFunctor>::
    warpImage(Image<InputFormat> const& inputImage,
              typename Image<OutputFormat>::PixelType defaultValue) const
    {
      if((inputImage.rows() != m_inputRows)
         || (inputImage.columns() != m_inputColumns)) {
        std::ostringstream message;
        message
          << "InputImage (" << inputImage.rows() << "x" << inputImage.columns()
          << ") doesn't match expected dimensions (" << m_inputRows << "x"
          << m_inputColumns << ").";
        BRICK_THROW(brick::common::ValueException, "ImageWarper::warpImage()",
                    message.str().c_str());
      }
      Image<OutputFormat> outputImage(
        m_lookupTable.rows(), m_lookupTable.columns());
      for(size_t ii = 0; ii < m_lookupTable.size(); ++ii) {
        SampleInfo const& sampleInfo = m_lookupTable(ii);
        if(sampleInfo.isInBounds) {
          typename Image<OutputFormat>::PixelType& outputPixel =
            outputImage[ii];
          size_t inputIndex = sampleInfo.index00;

          outputPixel = sampleInfo.c00 * inputImage[inputIndex];
          ++inputIndex;
          outputPixel += sampleInfo.c01 * inputImage[inputIndex];
          inputIndex += m_inputColumns;
          outputPixel += sampleInfo.c11 * inputImage[inputIndex];
          --inputIndex;
          outputPixel += sampleInfo.c10 * inputImage[inputIndex];
        } else {
          outputImage[ii] = defaultValue;
        }
      }
      return outputImage;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEWARPER_IMPL_HH */
