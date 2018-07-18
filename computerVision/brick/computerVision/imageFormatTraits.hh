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
      static PixelType getZeroPixel() {return PixelType();}
      static bool isIntegral() {return true;}
    };


    /**
     ** This class lets you map from pixel component type and number
     ** of channels back to the relevant image format.  For example,
     ** if you needed a single channel image with pixels that are
     ** UnsignedInt16, you could use
     ** ImageFormatTraitsIdentifierGray<UnsignedInt16>::Format to see
     ** that the relevant ImageFormat is GRAY16.
     **/
    template <class Type>
    class ImageFormatIdentifierGray {
    public:
      // Default case is not implemented.  Specializations are in
      // imageFormatTraits_impl.hh.
      static ImageFormat const Format = NO_FORMAT;
    };

  } // namespace computerVision

} // namespace brick


// Include file containing template specializations.
#include <brick/computerVision/imageFormatTraits_impl.hh>

#endif /* ifndef BRICK_COMPUTERVISION_IMAGEFORMATTRAITS_HH */
