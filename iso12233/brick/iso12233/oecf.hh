/**
***************************************************************************
* @file brick/iso12233/oecf.hh
*
* Header file declaring example functors that can be passed as the OECF
* argument of brick::iso12233::iso12233().
*
* Copyright (C) 2020 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_ISO12233_OECF_HH
#define BRICK_ISO12233_OECF_HH

#include <brick/iso12233/symbolImports.hh>

namespace brick {

  namespace iso12233 {

    /**
     * All this functor does is shift the pixel values in the image
     * patch so that the minimum value is 0.  This removes any extra
     * DC component from the patch, and prevents this DC component
     * from artificially depressing the MTF computed by iso12233().
     */
    template <class FloatType, ImageFormat Format>
    class ShiftOECF {
    public:

      typedef typename
      brick::computerVision::ImageFormatTraits<Format>::PixelType
      PixelType;
      
      /**
       * Constructor.
       *
       * @param inputPatch The image patch to be normalized.
       */
      ShiftOECF(Image<Format> const& inputPatch);
      
      
      /**
       * Destructor.
       */
      ~ShiftOECF() {}


      /**
       * Application operator.
       *
       * @param pixel The pixel value to be normalized.
       * @return The return value is the normalize pixel value.
       */
      inline FloatType
      operator()(PixelType const& pixel) const;

      static inline FloatType
      convertToFloat(PixelType const& pixel);

    private:
      FloatType m_minimumValue;
    };
    

  } // namespace iso12233

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/iso12233/oecf_impl.hh>

#endif /* #ifndef BRICK_ISO12233_OECF_HH */
