/**
***************************************************************************
* @file brick/computerVision/imagePyramidBinomial.hh
*
* Header file declaring a class template for constructing scale-space
* image pyramids.
*
* Copyright (C) 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEPYRAMIDBINOMIAL_HH
#define BRICK_COMPUTERVISION_IMAGEPYRAMIDBINOMIAL_HH

#include <stdint.h>
#include <deque>

#include <brick/computerVision/image.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This class template generates a scale-space image pyramid by
     ** repeatedly subsampling the input image.  Unlike the
     ** imagePyramid class template, this pyramid is constrained to
     ** use one image per octave, and the image filtering method
     ** (along each axis) is limited to a 3 element binomial kernel.
     **
     ** Template argument Format specifies the format of the input
     ** image.  InternalFormat specifies the image type used to
     ** accumulate results during low-pass filtering.
     **
     ** Example usage:
     **
     ** @code
     **   Image<GRAY8> inputImage = readPGM8(inputFileName);
     **   ImagePyramidBinomial<GRAY8, GRAY16> pyramid(inputImage);
     **   Image<GRAY8> smallestImage = pyramid.getLevel(
     **     pyramid.getNumberOfLevels() - 1);
     ** @endcode
     **/
    template <ImageFormat Format, ImageFormat InternalFormat>
    class ImagePyramidBinomial {
    public:

      // ========= Public member functions. =========

      /**
       * Construct an image pyramid from an input image.  The input
       * image will be stored as the base of the pyramid, then
       * low-pass filtered and to create subsequent levels.  If
       * desired, low-pass filtered images will be subtracted from
       * unfiltered images to create approximation of the
       * Laplacian-of-Gaussian scale space at each level except the
       * last.  The last level of the pyramid will simply be low-pass
       * filtered, regardless of whether or not the
       * Laplacian-of-Gaussian approximation is enabled.
       *
       * WARNING: the filtering currently zeros out the edges of the
       * low-pass-filtered image (where the filter doesn't fit on the
       * image). This creates a high-frequency edge several pixels
       * from the image borders.  This also makes a border of
       * un-subtracted pixels in the Laplacian-of-Gaussian
       * approximation.  The width of this border is reported by
       * member function getBorderSize().
       *
       * @param inputImage This argument is the image to be
       * downsampled.  It is either shallow-copied or deep-copied into
       * the base of the pyramid, depending on the value of
       * constructor argument isDeepCopyImage.
       *
       * @param levels This argument specifies how many pyramid levels
       * should be created.  Sitting this argument to zero tells the
       * pyramid to generate as many levels as possible.  Regardless
       * of whether or not levels is set to zero, the actual number of
       * pyramid levels generated is, subject to the limit imposed by
       * argument minimumImageSize.  For example, if levels is set to
       * 10, but only 3 pyramid levels can be generated without
       * resulting in an image size smaller than minimumImageSize,
       * then only 3 pyramid levels will be generated.
       *
       * @param minimumImageSize This argument specifies a lower limit
       * to the size of the smallest image of the pyramid.  Levels
       * will be added to the pyramid only to the extent that the
       * smallest image is larger than this in both axes.
       *
       * @param isBandPass Setting this argument to true enables the
       * Difference-of-Gaussians approximation to Laplacian of
       * Guassian scale space.  If this argument is false, only
       * low-pass-filtered (and then subsampled) pyramid levels will
       * be computed.
       *
       * @param isDeepCopyImage If this parameter is set to true, then
       * the image will be deep copied into the base of the image
       * pyramid.  Otherwise, the image will be shallow-copied (with
       * reference count), and the base of the pyramid will use the
       * same memory as the input image.
       */
      ImagePyramidBinomial(Image<Format> const& inputImage,
                           uint32_t levels = 0,
                           uint32_t minimumImageSize = 6,
                           bool isBandPass = true,
                           bool isDeepCopyImage = true);


      /**
       * This member function accepts pixel coordinates at one level
       * of the image and returns the matching coordinates in any
       * other level.  Pixel [0, 0] is in the upper left corner of the
       * image.  When converting from an earlier (larger image) level
       * to a later level, it is possible to input coordinates that
       * fall between pixels on the smaller image.  In this case, the
       * resulting coordinates will be rounded down.
       *
       * @param inCoords This argument specifies the coordinates to be
       * converted.
       *
       * @param fromLevel This argument specifies what pyramid level
       * the input coordinates apply to.
       *
       * @param toLevel This argument specifies to which pyramid level
       * the input coordinates should be converted.
       *
       * @return The return value is the converted coordinate pair.
       */
      brick::numeric::Index2D
      convertImageCoordinates(brick::numeric::Index2D const& inCoords,
                              unsigned int fromLevel,
                              unsigned int toLevel);


      /**
       * This member function accepts pixel coordinates at one level
       * of the image and returns the matching coordinates in any
       * other level.  Pixel coordinate (0, 0) is the upper left
       * corner of the upper left pixel (at each level).  The center
       * of the upper left pixel is (0.5, 0.5), and the X coordinate
       * increases along the top row of the image.
       *
       * @param inCoords This argument specifies the coordinates to be
       * converted.
       *
       * @param fromLevel This argument specifies what pyramid level
       * the input coordinates apply to.
       *
       * @param toLevel This argument specifies to which pyramid level
       * the input coordinates should be converted.
       *
       * @return The return value is the converted coordinate pair.
       */
      template <class Type>
      brick::numeric::Vector2D<Type>
      convertImageCoordinates(brick::numeric::Vector2D<Type> const& inCoords,
                              unsigned int fromLevel,
                              unsigned int toLevel);


      /**
       * Returns the number of columns at each of the left and right
       * borders of the image in which the filtered image data is
       * invalid due to the filter kernel overlapping the edge of the
       * image.
       *
       * @return The return value is the width, in pixels, of the
       * invalid region.
       */
      unsigned int
      getBorderSizeLeftRight();


      /**
       * Returns the number of rows at each of the top and bottom
       * borders of the image in which the filtered image data is
       * invalid due to the filter kernel overlapping the edge of the
       * image.
       *
       * @return The return value is the heigh, in pixels, of the
       * invalid region.
       */
      unsigned int
      getBorderSizeTopBottom();


      /**
       * This member function returns a shallow-copy of the specified
       * pyramid level.  If you need a deep copy, you can use
       * "myPyramid.getLevel().copy()."  The returned image is a
       * low-pass filtered and subsampled version of the original.
       *
       * @param levelIndex This argument is the requested level.
       * Level 0 is full resolution, with subsequent levels being more
       * and more subsampled.
       *
       * @return The return value is a shallow-copy of the requested
       * image.
       */
      Image<Format>&
      getLevel(unsigned int levelIndex);


      /**
       * Returns the number of levels in the pyramid.
       *
       * @return The return value specifies how many levels there are
       * in the pyramid.
       */
      unsigned int
      getNumberOfLevels();


    private:

      typedef typename ImageFormatTraits<Format>::PixelType
        PixelType;
      typedef typename ImageFormatTraits<InternalFormat>::PixelType
        InternalPixelType;


      void
      filterColumns(
        brick::numeric::Array1D<PixelType> outputRow,
        std::deque< brick::numeric::Array1D<InternalPixelType> > const&
          inputRowBuffer);

      Image<Format>
      filterImage(Image<Format> const& inputImage);

      brick::numeric::Array1D<InternalPixelType>
      filterRow(brick::numeric::Array1D<PixelType> const& inputRow);

      Image<Format>
      filterAndSubsampleImage(Image<Format> const& inputImage);

      brick::numeric::Array1D<InternalPixelType>
      filterAndSubsampleRow(brick::numeric::Array1D<PixelType> const& inputRow);

      Image<Format>
      subsampleImage(Image<Format> const& inputImage);

      int m_borderSizeLeftRight;
      int m_borderSizeTopBottom;
      std::vector< Image<Format> > m_pyramid;
    };

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/imagePyramidBinomial_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEPYRAMIDBINOMIAL_HH */
