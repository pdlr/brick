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

      ImagePyramid(Image<Format> const& inputImage,
                   double scaleFactorPerLevel = 2.0,
                   unsigned int levels = 0,
                   bool isBandPass = true);


      brick::numeric::Vector2D<double>
      convertImageCoordinates(brick::numeric::Vector2D<double> const& inCoords,
                              unsigned int fromLevel,
                              unsigned int toLevel);

      
      Image<Format>&
      getLevel(unsigned int levelIndex);


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
