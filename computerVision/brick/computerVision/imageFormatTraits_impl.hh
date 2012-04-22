/**
***************************************************************************
* @file brick/computerVision/imageFormatTraits_impl.hh
*
* Header file specializing the ImageFormatTraits class template.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEFORMATTRAITS_IMPL_HH
#define BRICK_COMPUTERVISION_IMAGEFORMATTRAITS_IMPL_HH

// This file is included by imageFormatTraits.hh, and should not be
// directly included by user code, so no need to include
// imageFormatTraits.hh here.
// 
// #include <brick/computerVision/imageFormatTraits.hh>

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
  
    template<>
    class ImageFormatTraits<GRAY1> {
    public:
      typedef bool PixelType;
      typedef bool ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return true;}
    };
    

    template<>
    class ImageFormatTraits<GRAY8> {
    public:
      typedef brick::common::UnsignedInt8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return true;}
    };
    

    template<>
    class ImageFormatTraits<GRAY16> {
    public:
      typedef brick::common::UnsignedInt16 PixelType;
      typedef brick::common::UnsignedInt16 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return true;}
    };
    

    template<>
    class ImageFormatTraits<GRAY32> {
    public:
      typedef brick::common::UnsignedInt32 PixelType;
      typedef brick::common::UnsignedInt32 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return true;}
    };
    

    template<>
    class ImageFormatTraits<GRAY64> {
    public:
      typedef brick::common::UnsignedInt64 PixelType;
      typedef brick::common::UnsignedInt64 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return true;}
    };
    

    template<>
    class ImageFormatTraits<GRAY_SIGNED16> {
    public:
      typedef brick::common::Int16 PixelType;
      typedef brick::common::Int16 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return true;}
    };
    

    template<>
    class ImageFormatTraits<GRAY_SIGNED32> {
    public:
      typedef brick::common::Int32 PixelType;
      typedef brick::common::Int32 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return true;}
    };
    

    template<>
    class ImageFormatTraits<GRAY_FLOAT32> {
    public:
      typedef brick::common::Float32 PixelType;
      typedef brick::common::Float32 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return false;}
    };
    

    template<>
    class ImageFormatTraits<GRAY_FLOAT64> {
    public:
      typedef brick::common::Float64 PixelType;
      typedef brick::common::Float64 ComponentType;
      static size_t getNumberOfComponents() {return 1;}
      static bool isIntegral() {return false;}
    };
    

    template<>
    class ImageFormatTraits<RGB8> {
    public:
      typedef PixelRGB8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
      static bool isIntegral() {return true;}
    };
    

    template<>
    class ImageFormatTraits<RGB16> {
    public:
      typedef PixelRGB16 PixelType;
      typedef brick::common::UnsignedInt16 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
      static bool isIntegral() {return true;}
    };


    template<>
    class ImageFormatTraits<RGB_SIGNED16> {
    public:
      typedef PixelRGBSigned16 PixelType;
      typedef brick::common::Int16 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
      static bool isIntegral() {return true;}
    };

  
    template<>
    class ImageFormatTraits<RGB_SIGNED32> {
    public:
      typedef PixelRGBSigned32 PixelType;
      typedef brick::common::Int32 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
      static bool isIntegral() {return true;}
    };

  
    template<>
    class ImageFormatTraits<RGB_FLOAT32> {
    public:
      typedef PixelRGBFloat32 PixelType;
      typedef brick::common::Float32 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
      static bool isIntegral() {return false;}
    };


    template<>
    class ImageFormatTraits<RGB_FLOAT64> {
    public:
      typedef PixelRGBFloat64 PixelType;
      typedef brick::common::Float64 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
      static bool isIntegral() {return false;}
    };


    template<>
    class ImageFormatTraits<HSV_FLOAT64> {
    public:
      typedef PixelHSVFloat64 PixelType;
      typedef brick::common::Float64 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
      static bool isIntegral() {return false;}
    };
  
  
    template<>
    class ImageFormatTraits<YIQ_FLOAT64> {
    public:
      typedef PixelYIQFloat64 PixelType;
      typedef brick::common::Float64 ComponentType;
      static size_t getNumberOfComponents() {return 3;}
      static bool isIntegral() {return false;}
    };
  
  
    template<>
    class ImageFormatTraits<RGBA8> {
    public:
      typedef PixelRGBA8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 4;}
      static bool isIntegral() {return true;}
    };
    

    template<>
    class ImageFormatTraits<BGRA8> {
    public:
      typedef PixelBGRA8 PixelType;
      typedef brick::common::UnsignedInt8 ComponentType;
      static size_t getNumberOfComponents() {return 4;}
      static bool isIntegral() {return true;}
    };


    template <>
    class ImageFormatIdentifierGray<bool> {
    public:
      static ImageFormat const Format = GRAY1;
    };


    template <>
    class ImageFormatIdentifierGray<brick::common::UnsignedInt8> {
    public:
      static ImageFormat const Format = GRAY8;
    };


    template <>
    class ImageFormatIdentifierGray<brick::common::UnsignedInt16> {
    public:
      static ImageFormat const Format = GRAY16;
    };


    template <>
    class ImageFormatIdentifierGray<brick::common::UnsignedInt32> {
    public:
      static ImageFormat const Format = GRAY32;
    };


    template <>
    class ImageFormatIdentifierGray<brick::common::UnsignedInt64> {
    public:
      static ImageFormat const Format = GRAY64;
    };


    template <>
    class ImageFormatIdentifierGray<brick::common::Int16> {
    public:
      static ImageFormat const Format = GRAY_SIGNED16;
    };
    

    template <>
    class ImageFormatIdentifierGray<brick::common::Int32> {
    public:
      static ImageFormat const Format = GRAY_SIGNED32;
    };
    

    template <>
    class ImageFormatIdentifierGray<brick::common::Float32> {
    public:
      static ImageFormat const Format = GRAY_FLOAT32;
    };
    

    template <>
    class ImageFormatIdentifierGray<brick::common::Float64> {
    public:
      static ImageFormat const Format = GRAY_FLOAT64;
    };
    
  } // namespace computerVision

} // namespace brick

#endif /* ifndef BRICK_COMPUTERVISION_IMAGEFORMATTRAITS_IMPL_HH */
