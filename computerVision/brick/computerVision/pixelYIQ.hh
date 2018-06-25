/**
***************************************************************************
* @file brick/computerVision/pixelYIQ.hh
*
* Header file declaring the PixelYIQ class template.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_PIXELYIQ_HH
#define BRICK_COMPUTERVISION_PIXELYIQ_HH

#include <brick/common/types.hh>

namespace brick {

  namespace computerVision {

    template<class TYPE>
    struct PixelYIQ
    {
    public:

      /**
       * This constructor makes no guarantees about the color of the
       * pixel.
       */
      PixelYIQ()
        : luma(), inPhase(), quadrature() {}


      /**
       * This constructor explicitly sets the pixel value.
       *
       * @param lumaValue This argument specifies the luma value for
       * the new pixel.
       *
       * @param inPhaseValue This argument specifies the inPhase value
       * for the new pixel.
       *
       * @param quadratureValue This argument specifies the quadrature
       * value for the new pixel.
       */
      PixelYIQ(const TYPE& lumaValue, const TYPE& inPhaseValue,
               const TYPE& quadratureValue)
        : luma(lumaValue), inPhase(inPhaseValue), quadrature(quadratureValue) {}


      /**
       * The destructor deletes cleans up for deletion.
       */
      ~PixelYIQ() {}


      /**
       * This member function copies the pixel component values, in
       * order, from consecutive iterator targets, incrementing the
       * iterator after each copy.
       *
       * @param iter This argument is the iterator from which to copy the
       * pixel components.
       *
       * @return The return value is a reference to iter after the
       * copying and incrementing is done.
       */
      template<class Iter>
      inline Iter&
      copyFromIterator(Iter& iter);


      /**
       * This member function assigns the pixel component values, in
       * order, to consecutive iterator targets, incrementing the
       * iterator after each assignment.
       *
       * @param iter This argument is the iterator to which to copy the
       * pixel components.
       *
       * @return The return value is a reference to iter after the
       * copying and incrementing is done.
       */
      template<class Iter>
      inline Iter&
      copyToIterator(Iter& iter);


      /* ====== Public data members ====== */

      TYPE luma;
      TYPE inPhase;
      TYPE quadrature;


      /* ====== Static member functions ====== */

      /**
       * This static member function indicates whether or not a pixel
       * instance is memory identical to a contiguous array of Component
       * type.
       *
       * @return The return value is true if the pixel structure is not
       * padded by the compiler.
       */
      static inline bool
      isContiguous();

    };


    typedef PixelYIQ<brick::common::UnsignedInt8> PixelYIQ8;
    typedef PixelYIQ<brick::common::UnsignedInt16> PixelYIQ16;
    typedef PixelYIQ<brick::common::Int16> PixelYIQSigned16;
    typedef PixelYIQ<brick::common::Int32> PixelYIQSigned32;
    typedef PixelYIQ<brick::common::Float32> PixelYIQFloat32;
    typedef PixelYIQ<brick::common::Float64> PixelYIQFloat64;


    /**
     * This operator subtracts the values of the individual color
     * components of its arguments.
     *
     * @param pixel0 The color component values of pixel1 will be
     * subtracted from the color component values of this pixel.
     *
     * @param pixel1 The color component values of this pixel will be
     * subtracted from the color component values of pixel1.
     *
     * @return The return value is a pixel in which each color component
     * value is the difference of the corresponding values in the two
     * input pixels.
     */
    template<class TYPE>
    inline PixelYIQ<TYPE>
    operator-(const PixelYIQ<TYPE>& pixel0, const PixelYIQ<TYPE>& pixel1);


    /**
     * This operator returns true if the contents of the two argments
     * are identical, false otherwise.
     *
     * @param pixel0 This argument is the first pixel value to be compared.
     *
     * @param pixel1 This argument is the second pixel value to be compared.
     *
     * @return The return value indicates whether the two pixels have
     * identical values.
     */
    template<class TYPE>
    inline bool
    operator==(const PixelYIQ<TYPE>& pixel0, const PixelYIQ<TYPE>& pixel1);

  } // namespace computerVision

} // namespace brick

/* ============ Definitions of inline & template functions ============ */


namespace brick {

  namespace computerVision {

    // This member function copies the pixel component values, in
    // order, from consecutive iterator targets, incrementing the
    // iterator after each copy.
    template<class TYPE>
    template<class Iter>
    inline Iter&
    PixelYIQ<TYPE>::
    copyFromIterator(Iter& iter)
    {
      luma = *(iter++);
      inPhase = *(iter++);
      quadrature = *(iter++);
      return iter;
    }


    // This member function assigns the pixel component values, in
    // order, to consecutive iterator targets, incrementing the
    // iterator after each assignment.
    template<class TYPE>
    template<class Iter>
    inline Iter&
    PixelYIQ<TYPE>::
    copyToIterator(Iter& iter)
    {
      *(iter++) = luma;
      *(iter++) = inPhase;
      *(iter++) = quadrature;
      return iter;
    }


    // This static member function indicates whether or not a pixel
    // instance is memory identical to a contiguous array of Component
    // type.
    template<class TYPE>
    inline bool
    PixelYIQ<TYPE>::
    isContiguous()
    {
      TYPE dummy;
      TYPE* typePtr = &dummy;
      PixelYIQ<TYPE>* pixelPtr = reinterpret_cast<PixelYIQ<TYPE>*>(typePtr);
      PixelYIQ<TYPE>* pixelPtr2 = pixelPtr + 1;
      return
        ((reinterpret_cast<TYPE*>(&(pixelPtr2->luma)) == &(typePtr[3]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->inPhase)) == &(typePtr[4]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->quadrature))
             == &(typePtr[5])));
    }


    // This operator subtracts the values of the individual color
    // components of its arguments.
    template<class TYPE>
    inline PixelYIQ<TYPE>
    operator-(const PixelYIQ<TYPE>& pixel0, const PixelYIQ<TYPE>& pixel1)
    {
      return PixelYIQ<TYPE>(
        pixel0.luma - pixel1.luma,
        pixel0.inPhase - pixel1.inPhase,
        pixel0.quadrature - pixel1.quadrature);
    }


    // This operator returns true if the contents of the two argments
    // are identical, false otherwise.
    template<class TYPE>
    inline bool
    operator==(const PixelYIQ<TYPE>& pixel0, const PixelYIQ<TYPE>& pixel1)
    {
      return (pixel0.luma == pixel1.luma
              && pixel0.inPhase == pixel1.inPhase
              && pixel0.quadrature == pixel1.quadrature);
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_PIXELYIQ_HH */
