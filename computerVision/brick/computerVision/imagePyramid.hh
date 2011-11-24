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

#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This class template generates a scale-space image pyramid by
     ** repeatedly subsampling the input image.
     **
     ** Template argument Format specifies the format of the input
     ** image.
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
       * unfiltered images to create a Laplacian-of-Gaussian
       * approximation at each level except the last.
       *
       * NOTE: the default filter kernel is quite slow.  If
       * performance is at all a concern, consider supplying your own.
       *
       * WARNING: the filtering currently zeros out the edges of the
       * low-pass-filtered image (where the filter doesn't fit on the
       * image). This creates a high-frequency edge several pixels
       * from the image borders.  This also makes a border of
       * un-subtracted pixels in the Laplacian of Gaussian images.
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
       * Laplacian of Guassian approximation.
       */
      ImagePyramid(Image<Format> const& inputImage,
                   double scaleFactorPerLevel = 2.0,
                   unsigned int levels = 0,
                   bool isBandPass = true);


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
       * This member function returns a shallow-copy of the specified
       * pyramid level.  If you need a deep copy, you can use
       * "myPyramid.getLevel().copy()."
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

      bool
      isIntegral(double scaleFactor, int& integerScaleFactor);


      Image<Format>
      subsampleImage(Image<Format> const& inputImage,
                     double scaleFactor);


      std::vector< Image<Format> > m_pyramid;
    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/imagePyramid_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEPYRAMID_HH */
