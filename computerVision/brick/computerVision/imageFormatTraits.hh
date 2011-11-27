/**
***************************************************************************
* @file brick/computerVision/imageFormatTraits.hh
*
* Header file declaring ImageFormatTraits classes.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEFORMATTRAITS_HH
#define BRICK_COMPUTERVISION_IMAGEFORMATTRAITS_HH

#include <stdlib.h>
#include <brick/common/types.hh>
#include <brick/computerVision/imageFormat.hh>
#include <brick/computerVision/pixelBGRA.hh>
#include <brick/computerVision/pixelHSV.hh>
#include <brick/computerVision/pixelRGB.hh>
#include <brick/computerVision/pixelRGBA.hh>
#include <brick/computerVision/pixelYIQ.hh>

namespace brick {

  namespace computerVision {
  
    /**
     ** The ImageFormatTraits class template specifies the
     ** characteristics of the available image formats.
     ** Specializations for each value of ImageFormat should be placed
     ** in imageFormatTraits_impl.hh.
     **/
    template <ImageFormat>
    class ImageFormatTraits {
    public:
      typedef brick::common::UnsignedInt8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return true;}
    };

  } // namespace computerVision

} // namespace brick


// Include file containing template specializations.
#include <brick/computerVision/imageFormatTraits_impl.hh>

#endif /* ifndef BRICK_COMPUTERVISION_IMAGEFORMATTRAITS_HH */
