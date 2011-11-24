/**
***************************************************************************
* @file brick/computerVision/colorspaceConverter.hh
*
* Header file declaring ColorspaceConverter class template.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_COLORSPACECONVERTER_HH
#define BRICK_COMPUTERVISION_COLORSPACECONVERTER_HH

#include <functional>
#include <brick/computerVision/imageFormat.hh>
#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {

    template<ImageFormat FORMAT0, ImageFormat FORMAT1>
    class ColorspaceConverter
      : public std::unary_function<typename Image<FORMAT0>::PixelType,
                                   typename Image<FORMAT1>::PixelType>
    {
    public:
    
      /** 
       * The default constructor simply dispatches to the
       * std::unary_function constructor.
       */
      ColorspaceConverter()
        : std::unary_function<typename Image<FORMAT0>::PixelType,
                              typename Image<FORMAT1>::PixelType>()
        {}


      /** 
       * The destructor destroys the class instnace and cleans up any
       * resources.
       * 
       * @return The return value 
       */
      virtual ~ColorspaceConverter()
        {}


      /** 
       * The application operator does the format conversion for one
       * pixel.  We somewhat awkwardly dispatch to the two-argument
       * application operator because of the unfortunate way we've had
       * to declare the non-scalar pixel types.  The non-scalar pixel
       * types are declared as structs with C linkage so we can
       * guarantee their layout in memory.  I don't know of a way to
       * return a filled-in C struct which takes advantage of the return
       * value optimization, so we provide the two-argument version of
       * operator()() in order to sidestep the problem.  The
       * single-argument version dispatches to the two-argument version
       * to avoid code duplication, which makes it slightly slower than
       * it could be for scalar types.
       *
       * @param inputPixel The pixel to be converted.
       * 
       * @return The return value is the converted pixel.
       */
      inline
      typename Image<FORMAT1>::PixelType
      operator()(const typename Image<FORMAT0>::PixelType& inputPixel) {
        typename Image<FORMAT1>::PixelType outputPixel;
        this->operator()(inputPixel, outputPixel);
        return outputPixel;
      }
    

      /** 
       * The application operator does the format conversion for one
       * pixel, returning the result through a reference argument.  This
       * is the fastest way to do a pixel conversion.
       * 
       * @param inputPixel The pixel to be converted.
       * 
       * @param outputPixel This reference argument is set to the result
       * of the conversion.
       */
      inline
      void
      operator()(const typename Image<FORMAT0>::PixelType& inputPixel,
                 typename Image<FORMAT1>::PixelType& outputPixel) {
        // Default rule should work for many format combinations.
        // We'll specialize for the rest in the implementation file.
        // If you try to convert between formats for which this
        // static_cast isn't appropriate, and for which we haven't
        // written a specialization, then you'll either get a compile
        // error or an unexpected conversion result.
        outputPixel =
          static_cast<typename Image<FORMAT1>::PixelType>(inputPixel);
      }
    
    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/colorspaceConverter_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_COLORSPACECONVERTER_HH */
