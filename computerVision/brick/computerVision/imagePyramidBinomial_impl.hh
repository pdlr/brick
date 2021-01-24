/**
***************************************************************************
* @file brick/computerVision/imagePyramidBinomial_impl.hh
*
* Header file defining a class template for constructing scale-space
* image pyramids.
*
* Copyright (C) 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEPYRAMIDBINOMIAL_IMPL_HH
#define BRICK_COMPUTERVISION_IMAGEPYRAMIDBINOMIAL_IMPL_HH

// This file is included by imagePyramidBinomial.hh, and should not be
// directly included by user code, so no need to include
// imagePyramidBinomial.hh here.
//
// #include <brick/computerVision/imagePyramidBinomial.hh>

#include <brick/computerVision/pixelOperations.hh>
#include <brick/computerVision/utilities.hh>
#include <brick/numeric/numericTraits.hh>

namespace brick {

  namespace computerVision {

    template <ImageFormat Format, ImageFormat InternalFormat>
    ImagePyramidBinomial<Format, InternalFormat>::
    ImagePyramidBinomial(Image<Format> const& inputImage,
                         uint32_t levels,
                         uint32_t minimumImageSize,
                         bool isBandPass,
                         bool isDeepCopyImage)
      : m_borderSizeLeftRight(-1),
        m_borderSizeTopBottom(-1),
        m_pyramid()
    {
      // How many levels are implied by minimumImageSize?
      uint32_t automaticLevels = 0;
      {
        // How many pyramid levels?  Well, enough that we get close to
        // minimumImageSize, but no smaller.  This implies that
        // minimumImageSize * 2^numberOfLevels is less than or equal
        // to input image size, but minimumImageSize *
        // scaleFactorPerLevel^(numberOfLevels + 1) is greater than
        // input image size.  That is, numberOfLevels is the largest
        // integer less or equal to the variable nn in the following equation:
        //
        //   minimumImageSize * 2^nn = inImageSize
        //
        //   2^nn = inImageSize / minimumImageSize
        //
        //   nn = ln(inImageSize / minimumImageSize) / ln(2)
        uint32_t limitingDimension =
          std::min(inputImage.rows(), inputImage.columns());
        double targetSizeRatio = (static_cast<double>(limitingDimension)
                                  / static_cast<double>(minimumImageSize));
        automaticLevels = std::max(
          static_cast<uint32_t>(std::log(targetSizeRatio) / std::log(2.0)),
          automaticLevels);
      }

      // If user didn't supply a non-zero value for levels, we're done.
      if(levels == 0) {
        levels = automaticLevels;
      }

      // Mediate between the user specified and automatically
      // generated pyramid sizes.
      levels = std::min(levels, automaticLevels);

      // The next two lines rely on the fact that our binomial kernel
      // is 3x3.
      m_borderSizeLeftRight = 1;
      m_borderSizeTopBottom = 1;

      // Start off the pyramid.
      Image<Format> currentImage;
      if(isDeepCopyImage) {
        currentImage = inputImage.copy();
      } else {
        currentImage = inputImage;
      }
      m_pyramid.push_back(currentImage);
      --levels;

      while(0 != levels) {
        if(isBandPass) {
          Image<Format> filteredImage = this->filterImage(currentImage);
          m_pyramid[m_pyramid.size() - 1] -= filteredImage;
          currentImage = this->subsampleImage(filteredImage);
        } else {
          currentImage = this->filterAndSubsampleImage(currentImage);
        }
        m_pyramid.push_back(currentImage);
        --levels;
      }
    }


    // This member function accepts pixel coordinates at one level
    // of the image and returns the matching coordinates in any
    // other level.
    template <ImageFormat Format, ImageFormat InternalFormat>
    brick::numeric::Index2D
    ImagePyramidBinomial<Format, InternalFormat>::
    convertImageCoordinates(brick::numeric::Index2D const& inCoords,
                            unsigned int fromLevel,
                            unsigned int toLevel)
    {
      int row = inCoords.getRow();
      int column = inCoords.getColumn();
      if(fromLevel > toLevel) {
        row << fromLevel - toLevel;
        column << fromLevel - toLevel;
      } else if(fromLevel < toLevel) {
        row >> toLevel - fromLevel;
        column >> toLevel - fromLevel;
      }
      return brick::numeric::Index2D(row, column);
    }


    // This member function accepts pixel coordinates at one level
    // of the image and returns the matching coordinates in any
    // other level.
    template <ImageFormat Format, ImageFormat InternalFormat>
    template <class Type>
    brick::numeric::Vector2D<Type>
    ImagePyramidBinomial<Format, InternalFormat>::
    convertImageCoordinates(brick::numeric::Vector2D<Type> const& inCoords,
                            unsigned int fromLevel,
                            unsigned int toLevel)
    {
      Type row = inCoords.y();
      Type column = inCoords.x();
      while(fromLevel > toLevel) {
        row *= Type(2.0);
        column *= Type(2.0);
        --fromLevel;
      }
      while(fromLevel < toLevel) {
        row /= Type(2.0);
        column /= Type(2.0);
        ++fromLevel;
      }
      return brick::numeric::Vector2D<Type>(row, column);
    }


    template <ImageFormat Format, ImageFormat InternalFormat>
    unsigned int
    ImagePyramidBinomial<Format, InternalFormat>::
    getBorderSizeLeftRight()
    {
      if(m_borderSizeLeftRight < 0) {
        BRICK_THROW(brick::common::StateException,
                    "ImagePyramidBinomial::getBorderSizeLeftRight()",
                    "Border size has not been set yet.");
      }
      return static_cast<unsigned int>(m_borderSizeLeftRight);
    }


    template <ImageFormat Format, ImageFormat InternalFormat>
    unsigned int
    ImagePyramidBinomial<Format, InternalFormat>::
    getBorderSizeTopBottom()
    {
      if(m_borderSizeTopBottom < 0) {
        BRICK_THROW(brick::common::StateException,
                    "ImagePyramidBinomial::getBorderSizeLeftRight()",
                    "Border size has not been set yet.");
      }
      return static_cast<unsigned int>(m_borderSizeTopBottom);
    }


    template <ImageFormat Format, ImageFormat InternalFormat>
    Image<Format>&
    ImagePyramidBinomial<Format, InternalFormat>::
    getLevel(unsigned int levelIndex)
    {
      if(levelIndex >= m_pyramid.size()) {
        std::ostringstream message;
        message << "Argument levelIndex (" << levelIndex << ") "
                << "is invalid for a " << m_pyramid.size()
                << " level image pyramid.";
        BRICK_THROW(brick::common::IndexException,
                    "ImagePyramidBinomial::getLevel()", message.str().c_str());
      }
      return m_pyramid[levelIndex];
    }


    template <ImageFormat Format, ImageFormat InternalFormat>
    unsigned int
    ImagePyramidBinomial<Format, InternalFormat>::
    getNumberOfLevels()
    {
      return m_pyramid.size();
    }


    // ============== Private member functions below this line ==============

    template <ImageFormat Format, ImageFormat InternalFormat>
    void
    ImagePyramidBinomial<Format, InternalFormat>::
    filterColumns(
      brick::numeric::Array1D<
        typename ImagePyramidBinomial<Format, InternalFormat>::PixelType>
        outputRow,
      std::deque< brick::numeric::Array1D<
        typename ImagePyramidBinomial<Format, InternalFormat>::InternalPixelType> >
        const& inputRowBuffer)
    {
      if(inputRowBuffer.size() != 3) {
        BRICK_THROW(brick::common::LogicException,
                    "ImagePyramidBinomial::filterColumns()",
                    "Input row deque must have exaxtly 3 elements.");
      }
      if(outputRow.size() != inputRowBuffer[0].size()) {
        outputRow.reinit(inputRowBuffer[0].size());
      }

      for(uint32_t ii = 0; ii < outputRow.size(); ++ii) {
        InternalPixelType tempPixel =
          inputRowBuffer[0][ii]
          + multiplyPixel<2>(inputRowBuffer[1][ii])
          + inputRowBuffer[2][ii];
        InternalPixelType outputPixelValue = dividePixel<4>(tempPixel);
        outputRow[ii] = static_cast<PixelType>(outputPixelValue);
      }
    }


    template <ImageFormat Format, ImageFormat InternalFormat>
    Image<Format>
    ImagePyramidBinomial<Format, InternalFormat>::
    filterImage(Image<Format> const& inputImage)
    {
      // Check arguments.
      if(inputImage.rows() < 3 || inputImage.columns() < 3) {
        BRICK_THROW(brick::common::LogicException,
                    "ImagePyramidBinomial::filterImage()",
                    "Input image size is smaller than 3x3.  There must be "
                    "some unchecked arguments upstream.");
      }

      // Create an output image of the appropriate size.
      Image<Format> outputImage(inputImage.rows(), inputImage.columns());

      // Zero the first output row (where there will be no filtered data).
      std::fill(outputImage.getRow(0).begin(), outputImage.getRow(0).end(),
                PixelType(0));

      // Now iterate through all input rows, low pass filtering.
      uint32_t inputRowIndex = 0;
      uint32_t outputRowIndex = 0;
      std::deque< brick::numeric::Array1D<InternalPixelType> > rowBuffer;

      // First few rows just get low pass filtered, and then queued up
      // until we have enough rows in the pipeline to start generating
      // output.
      while(inputRowIndex < 2) {
        rowBuffer.push_back(this->filterRow(inputImage.getRow(inputRowIndex)));
        ++inputRowIndex;
      }

      // Generate output until we reach the bottom of the input image.
      while(inputRowIndex < inputImage.rows()) {
        rowBuffer.push_back(this->filterRow(inputImage.getRow(inputRowIndex)));
        this->filterColumns(outputImage.getRow(outputRowIndex), rowBuffer);
        rowBuffer.pop_front();
        ++outputRowIndex;
        ++inputRowIndex;
      }

      // Zero the last output row (where there is no filtered data).
      while(outputRowIndex < outputImage.rows()) {
        std::fill(outputImage.getRow(outputRowIndex).begin(),
                  outputImage.getRow(outputRowIndex).end(),
                  PixelType(0));
        ++outputRowIndex;
      }

      return outputImage;
    }


    template <ImageFormat Format, ImageFormat InternalFormat>
    brick::numeric::Array1D<typename ImagePyramidBinomial<Format, InternalFormat>::InternalPixelType>
    ImagePyramidBinomial<Format, InternalFormat>::
    filterRow(brick::numeric::Array1D<PixelType> const& inputRow)
    {
      brick::numeric::Array1D<InternalPixelType> outputRow(inputRow.size());
      uint32_t sizeMinusOne = outputRow.size() - 1;

      outputRow[0] =
        ImageFormatTraits<InternalFormat>().getZeroPixel();

      for(uint32_t outputRowIndex = 1;
          outputRowIndex < sizeMinusOne;
          ++outputRowIndex) {
        InternalPixelType tempPixel =
          static_cast<InternalPixelType>(inputRow[outputRowIndex - 1])
          + multiplyPixel<2>(
            static_cast<InternalPixelType>(inputRow[outputRowIndex]))
          + static_cast<InternalPixelType>(inputRow[outputRowIndex + 1]);
        outputRow[outputRowIndex] = dividePixel<4>(tempPixel);

      }

      outputRow[sizeMinusOne] =
        ImageFormatTraits<InternalFormat>().getZeroPixel();

      return outputRow;
    }


    template <ImageFormat Format, ImageFormat InternalFormat>
    Image<Format>
    ImagePyramidBinomial<Format, InternalFormat>::
    filterAndSubsampleImage(Image<Format> const& inputImage)
    {
      // Check arguments.
      if(inputImage.rows() < 3 || inputImage.columns() < 3) {
        BRICK_THROW(brick::common::LogicException,
                    "ImagePyramidBinomial::filterAndSubsampleImage()",
                    "Input image size is smaller than 3x3.  There must be "
                    "some unchecked arguments upstream.");
      }

      // Create an output image of the appropriate size.
      Image<Format> outputImage(inputImage.rows() / 2,
                                inputImage.columns() / 2);

      // Zero the first output row (where there will be no filtered data).
      std::fill(outputImage.getRow(0).begin(), outputImage.getRow(0).end(),
                PixelType(0));

      // Now iterate through all input rows, low pass filtering and
      // subsampling..
      uint32_t inputRowIndex = 0;
      uint32_t outputRowIndex = 0;
      std::deque< brick::numeric::Array1D<InternalPixelType> > rowBuffer;

      // First row just gets low pass filtered & subsampled, and then
      // queued up until we have enough rows in the pipeline to start
      // generating output.
      rowBuffer.push_back(
        this->filterAndSubsampleRow(inputImage.getRow(inputRowIndex)));
      ++inputRowIndex;

      // Generate output until we reach the bottom of the input image.
      while(inputRowIndex < (inputImage.rows() - 1)) {
        // Remember we're subsampling, so go ahead and filter the next
        // input row without generating an output row.
        rowBuffer.push_back(
          this->filterAndSubsampleRow(inputImage.getRow(inputRowIndex)));
        ++inputRowIndex;

        // Now filter the next input row and generate an output row.
        rowBuffer.push_back(
          this->filterAndSubsampleRow(inputImage.getRow(inputRowIndex)));
        this->filterColumns(outputImage.getRow(outputRowIndex), rowBuffer);
        ++inputRowIndex;

        // On to the next row!
        rowBuffer.pop_front();
        rowBuffer.pop_front();
        ++outputRowIndex;
      }

      // If necessary, zero the last output row (where there is no
      // filtered data).
      while(outputRowIndex < outputImage.rows()) {
        std::fill(outputImage.getRow(outputRowIndex).begin(),
                  outputImage.getRow(outputRowIndex).end(),
                  PixelType(0));
        ++outputRowIndex;
      }

      return outputImage;
    }


    template <ImageFormat Format, ImageFormat InternalFormat>
    brick::numeric::Array1D<typename ImagePyramidBinomial<Format, InternalFormat>::InternalPixelType>
    ImagePyramidBinomial<Format, InternalFormat>::
    filterAndSubsampleRow(brick::numeric::Array1D<PixelType> const& inputRow)
    {
      brick::numeric::Array1D<InternalPixelType> outputRow(inputRow.size() / 2);
      uint32_t sizeMinusOne = inputRow.size() - 1;

      outputRow[0] =
        ImageFormatTraits<InternalFormat>().getZeroPixel();

      uint32_t outputRowIndex = 1;
      for(uint32_t inputRowIndex = 2;
          inputRowIndex < sizeMinusOne;
          inputRowIndex += 2) {
        InternalPixelType tempPixel =
          static_cast<InternalPixelType>(inputRow[inputRowIndex - 1])
          + multiplyPixel<2>(
              static_cast<InternalPixelType>(inputRow[inputRowIndex]))
          + static_cast<InternalPixelType>(inputRow[inputRowIndex + 1]);
        outputRow[outputRowIndex] = dividePixel<4>(tempPixel);
        ++outputRowIndex;
      }

      while(outputRowIndex < outputRow.size()) {
        outputRow[outputRowIndex] =
          ImageFormatTraits<InternalFormat>().getZeroPixel();
        ++outputRowIndex;
      }

      return outputRow;
    }


    template <ImageFormat Format, ImageFormat InternalFormat>
    Image<Format>
    ImagePyramidBinomial<Format, InternalFormat>::
    subsampleImage(Image<Format> const& inputImage)
    {
      Image<Format> outputImage(inputImage.rows() / 2,
                                inputImage.columns() / 2);
      uint32_t inputRowIndex = 0;
      for(uint32_t outputRowIndex = 0;
          outputRowIndex < outputImage.rows();
          ++outputRowIndex) {
        uint32_t inputColumnIndex = 0;
        for(uint32_t outputColumnIndex = 0;
            outputColumnIndex < outputImage.columns();
            ++outputColumnIndex) {
          outputImage(outputRowIndex, outputColumnIndex) =
            inputImage(inputRowIndex, inputColumnIndex);
          inputColumnIndex += 2;
        }
        inputRowIndex += 2;
      }
      return outputImage;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEPYRAMIDBINOMIAL_IMPL_HH */
