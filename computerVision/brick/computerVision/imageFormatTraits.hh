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
     ** characteristics of the available image formats.  Currently it
     ** specifies only the base type used to represent pixels.
     **/
    template <ImageFormat>
    class ImageFormatTraits {
    public:
      typedef brick::common::UnsignedInt8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };


    template<>
    class ImageFormatTraits<GRAY1> {
    public:
      typedef bool PixelType;
      typedef bool ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };
    

    template<>
    class ImageFormatTraits<GRAY8> {
    public:
      typedef brick::common::UnsignedInt8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };
    

    template<>
    class ImageFormatTraits<GRAY16> {
    public:
      typedef brick::common::UnsignedInt16 PixelType;
      typedef brick::common::UnsignedInt16 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };
    

    template<>
    class ImageFormatTraits<GRAY32> {
    public:
      typedef brick::common::UnsignedInt32 PixelType;
      typedef brick::common::UnsignedInt32 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };
    

    template<>
    class ImageFormatTraits<GRAY64> {
    public:
      typedef brick::common::UnsignedInt64 PixelType;
      typedef brick::common::UnsignedInt64 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };
    

    template<>
    class ImageFormatTraits<GRAY_SIGNED16> {
    public:
      typedef brick::common::Int16 PixelType;
      typedef brick::common::Int16 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };
    

    template<>
    class ImageFormatTraits<GRAY_SIGNED32> {
    public:
      typedef brick::common::Int32 PixelType;
      typedef brick::common::Int32 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };
    

    template<>
    class ImageFormatTraits<GRAY_FLOAT32> {
    public:
      typedef brick::common::Float32 PixelType;
      typedef brick::common::Float32 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };
    

    template<>
    class ImageFormatTraits<GRAY_FLOAT64> {
    public:
      typedef brick::common::Float64 PixelType;
      typedef brick::common::Float64 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
    };
    

    template<>
    class ImageFormatTraits<RGB8> {
    public:
      typedef PixelRGB8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
    };
    

    template<>
    class ImageFormatTraits<RGB16> {
    public:
      typedef PixelRGB16 PixelType;
      typedef brick::common::UnsignedInt16 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
    };


    template<>
    class ImageFormatTraits<RGB_SIGNED16> {
    public:
      typedef PixelRGBSigned16 PixelType;
      typedef brick::common::Int16 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
    };

  
    template<>
    class ImageFormatTraits<RGB_SIGNED32> {
    public:
      typedef PixelRGBSigned32 PixelType;
      typedef brick::common::Int32 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
    };

  
    template<>
    class ImageFormatTraits<RGB_FLOAT32> {
    public:
      typedef PixelRGBFloat32 PixelType;
      typedef brick::common::Float32 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
    };


    template<>
    class ImageFormatTraits<RGB_FLOAT64> {
    public:
      typedef PixelRGBFloat64 PixelType;
      typedef brick::common::Float64 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
    };


    template<>
    class ImageFormatTraits<HSV_FLOAT64> {
    public:
      typedef PixelHSVFloat64 PixelType;
      typedef brick::common::Float64 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
    };
  
  
    template<>
    class ImageFormatTraits<YIQ_FLOAT64> {
    public:
      typedef PixelYIQFloat64 PixelType;
      typedef brick::common::Float64 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
    };
  
  
    template<>
    class ImageFormatTraits<RGBA8> {
    public:
      typedef PixelRGBA8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 4;}
    };
    

    template<>
    class ImageFormatTraits<BGRA8> {
    public:
      typedef PixelBGRA8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 4;}
    };
    
  } // namespace computerVision

} // namespace brick

#endif /* ifndef BRICK_COMPUTERVISION_IMAGEFORMATTRAITS_HH */
