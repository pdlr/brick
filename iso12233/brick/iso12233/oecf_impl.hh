/**
***************************************************************************
* @file brick/iso12233/iso12233_impl.hh
*
* Header file defining functors that can be passed as the OECF
* argument of brick::iso12233::iso12233().
*
* Copyright (C) 2020 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_ISO12233_OECF_IMPL_HH
#define BRICK_ISO12233_OECF_IMPL_HH

// This file is included by oecf.hh, and should not be directly included
// by user code, so no need to include oecf.hh here.
//
// #include <brick/iso12233/oecf.hh>

namespace brick {

  namespace iso12233 {

    // Constructor.
    template <class FloatType, ImageFormat Format>
    ShiftOECF<FloatType, Format>::
    ShiftOECF(Image<Format> const& inputPatch)
      : m_minimumValue(std::numeric_limits<FloatType>::max())
    {
      for(auto pixel : inputPatch) {
        FloatType value = ShiftOECF<FloatType, Format>::convertToFloat(pixel);
        if(value < this->m_minimumValue) {
          this->m_minimumValue = value;
        }
      }
    }      
      

    // Application operator.
    template <class FloatType, ImageFormat Format>
    FloatType
    ShiftOECF<FloatType, Format>::
    operator()(PixelType const& pixel) const
    {
      return (ShiftOECF<FloatType, Format>::convertToFloat(pixel)
              - this->m_minimumValue);
    }


    // Default implementation should work for many image format /
    // float types.
    template <class FloatType, ImageFormat Format>
    FloatType
    ShiftOECF<FloatType, Format>::
    convertToFloat(PixelType const& pixel)
    {
      return static_cast<FloatType>(pixel);
    }
    
  } // namespace iso12233

} // namespace brick

#endif /* #ifndef BRICK_ISO12233_OECF_IMPL_HH */
