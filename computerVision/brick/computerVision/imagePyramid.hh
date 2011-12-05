/**
***************************************************************************
* @file brick/computerVision/imagePyramid.hh
*
* Header file declaring a class template for constructing scale-space
* image pyramids.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEPYRAMID_HH
#define BRICK_COMPUTERVISION_IMAGEPYRAMID_HH

#include <brick/computerVision/image.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This class template generates a scale-space image pyramid by
     ** repeatedly subsampling the input image.
     **
     ** Template argument Format specifies the format of the input
     ** image.  InternalFormat specifies the image type used to
     ** accumulate results during low-pass filtering, and KernelType
     ** specifies the element type of the filter kernel used to
     ** perform the low-pass filtering.
     **
     ** Example usage:
     **
     ** @code
     **   Image<GRAY8> inputImage = readPGM8(inputFileName);
     **   Image<GRAY_FLOAT32> floatImage = convertColorspace<GRAY_FLOAT32>(
     **     inputImage);
     **   ImagePyramid<GRAY_FLOAT32, GRAY_FLOAT32, common::Float32> pyramid(
     **     floatImage);
     **   Image<GRAY_FLOAT32> smallestImage = pyramid.getLevel(
     **     pyramid.getNumberOfLevels() - 1);
     ** @endcode
     **
     ** In this example, we choose a floating point format so that the
     ** result of the low-pass filtering operation is well behaved.
     ** Over time, filtering operations in brick::computerVision may
     ** mature to the point that discrete image formats (GRAY16,
     ** GRAY32) are a good (and more efficent) choice.
     **/
    template <ImageFormat Format, ImageFormat InternalFormat, class KernelType>
    class ImagePyramid {
    public:

      // ========= Public member functions. =========

      /** 
       * Construct an image pyramid using default filtering.  The
       * input image will be stored as the base of the pyramid, then
       * low-pass filtered and to create subsequent levels.  If
       * desired, low-pass filtered images will be subtracted from
       * unfiltered images to create a Difference-of-Gaussians
       * approximation to the Laplacian-of-Gaussian scale space at
       * each level except the last.  The last level of the pyramid
       * will simply be low-pass (Gaussian) filtered, regardless of
       * whether or not the Difference-of-Gaussians computation is
       * enabled.
       *
       * WARNING: the filtering currently zeros out the edges of the
       * low-pass-filtered image (where the filter doesn't fit on the
       * image). This creates a high-frequency edge several pixels
       * from the image borders.  This also makes a border of
       * un-subtracted pixels in the Difference-of-Gaussian images.
       * The width of this border is reported by member function
       * getBorderSize().
       * 
       * @param inputImage This argument is the image to be downsampled.
       * 
       * @param scaleFactorPerLevel This argument specifies the size
       * ratio between pyramid levels.  Integer scale factors are
       * detected and handled more efficiently than non-integer scale
       * factors.  The "special" scale factor 1.5 (special because
       * it's easy to do quickly) is not detected, but may be in a
       * subsequent release.
       * 
       * @param levels This argument specifies how many pyramid levels
       * should be created.  Sitting this argument to zero tells the
       * pyramid to generate as many levels as possible.
       * 
       * @param isBandPass Setting this argument to true enables the
       * Difference-of-Gaussians approximation to Laplacian of
       * Guassian scale space.  If this argument is false, only
       * low-pass-filtered (and then subsampled) pyramid levels will
       * be computed.
       */
      ImagePyramid(Image<Format> const& inputImage,
                   double scaleFactorPerLevel = 2.0,
                   unsigned int levels = 0,
                   bool isBandPass = true);


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
      brick::numeric::Vector2D<double>
      convertImageCoordinates(brick::numeric::Vector2D<double> const& inCoords,
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


      /** 
       * This member function returns a shallow-copy of the specified
       * scale-space pyramid level.  If you need a deep copy, you can
       * use "myPyramid.getScaleSpaceLevel().copy()."  If constructor
       * argument isBandPass was set to true, the returned image is a
       * band-pass (Difference of Gaussians) filtered and subsampled
       * version of the original.  Otherwise, the returned image is
       * empty.
       * 
       * @param levelIndex This argument is the requested level.
       * Level 0 is full resolution, with subsequent levels being more
       * and more subsampled.
       * 
       * @return The return value is a shallow-copy of the requested
       * image.
       */
      Image<Format>&
      getScaleSpaceLevel(unsigned int levelIndex);

    private:

      bool
      isIntegral(double scaleFactor, int& integerScaleFactor);


      Image<Format>
      subsampleImage(Image<Format> const& inputImage,
                     double scaleFactor);


      int m_borderSizeLeftRight;
      int m_borderSizeTopBottom;
      std::vector< Image<Format> > m_pyramid;
      double m_scaleFactorPerLevel;
    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/imagePyramid_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEPYRAMID_HH */
