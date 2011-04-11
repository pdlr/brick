/**
***************************************************************************
* @file brick/common/byteOrder.hh
*
* Header file declaring some useful routines related to endian-ness.
*
* (C) Copyright 2006-2007 David LaRose, dlr@alumni.carnegiemellon.edu
* See accompanying LICENSE file for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_BYTEORDER_HH
#define BRICK_COMMON_BYTEORDER_HH

#include <algorithm>
#include <sstream>
#include <brick/common/exception.hh>
#include <brick/common/types.hh>

namespace brick {

  namespace common {
    
    /**
     ** This enum provides a convenient way to represent the various
     ** machine-dependent byte orderings.
     **/
    enum ByteOrder {
      BRICK_BIG_ENDIAN,
      BRICK_LITTLE_ENDIAN,
    };


    /** 
     * This function returns the appropriate byte order for the platform
     * on which it is run.  For example, this function will return
     * BRICK_LITTLE_ENDIAN when run on a 386 machine.
     * 
     * @return The return value is the byte ordering used by the current
     * platform.
     */
    inline ByteOrder
    getByteOrder();


    /** 
     * This function swaps the byte order of a single value.
     * 
     * @param inputValue This argument is the value to be byte-swapped.
     * 
     * @param fromByteOrder This argument indicates the current byte
     * order of the values in the array.
     * 
     * @param toByteOrder This argument specifies the desired final byte
     * order for the data in the array.
     *
     * @return The byte-swapped value.
     */
    template <class Type>
    inline Type
    switchByteOrder(Type inputValue,
                    ByteOrder fromByteOrder,
                    ByteOrder toByteOrder);

    
    /** 
     * This function takes a pointer to a C-style array of values and
     * modifies the array in place so that it has a particular byte
     * order.
     * 
     * @param dataPtr This argument is a pointer to the C-style array of
     * values.
     * 
     * @param numberOfElements This argument indicates how many values
     * are in the array.
     * 
     * @param fromByteOrder This argument indicates the current byte
     * order of the values in the array.
     * 
     * @param toByteOrder This argument specifies the desired final byte
     * order for the data in the array.
     */
    template <class Type>
    inline void
    switchByteOrder(Type* dataPtr,
                    size_t numberOfElements,
                    ByteOrder fromByteOrder,
                    ByteOrder toByteOrder);

    
    /** 
     * This function takes a pointer to a C-style array of values and
     * copies it into another C-style array, swapping bytes if necessary
     * so that the output array has the specified byte order.
     * 
     * @param fromDataPtr This argument is a pointer to the C-style
     * array of input values.
     * 
     * @param numberOfElements This argument indicates how many values
     * are in both the input array and the output array.
     * 
     * @param toDataPtr This argument is a pointer to the C-style array
     * of output values.
     * 
     * @param fromByteOrder This argument indicates the current byte
     * order of the values in the array.
     * 
     * @param toByteOrder This argument specifies the desired final byte
     * order for the data in the array.
     */
    template <class Type>
    inline void
    switchByteOrder(const Type* fromDataPtr,
                    size_t numberOfElements,
                    Type* toDataPtr,
                    ByteOrder fromByteOrder,
                    ByteOrder toByteOrder);


  } // namespace common
  
} // namespace brick


/* ============ Implementation of template functions follows ============ */

#include <algorithm>

namespace brick {

  namespace common {

    /// @cond privateCode
    namespace privateCode {
    
      /**
       * This is an "under-the-hood" implementation function for the
       * switchByteOrder() function template.  It is not part of the
       * public interface, and may be removed later.
       *
       * The default implementation isn't smart enough to know how to switch
       * any byte ordering.  We use specializations for that.
       **/
      template <int Size>
      inline void
      genericSwitchByteOrder(UnsignedInt8* dataPtr,
                             size_t numberOfElements,
                             ByteOrder fromByteOrder,
                             ByteOrder toByteOrder)
      {
        if(fromByteOrder == toByteOrder) {
          return;
        }
        std::ostringstream message;
        message << "This function is not implemented for types of size"
                << Size << ".";
        BRICK_THROW(NotImplementedException, "genericSwitchByteOrder()",
                    message.str().c_str());
      }
  

      template <>
      inline void
      genericSwitchByteOrder<2>(UnsignedInt8* dataPtr,
                                size_t numberOfElements,
                                ByteOrder fromByteOrder,
                                ByteOrder toByteOrder)
      {
        if(fromByteOrder == toByteOrder) {
          return;
        }
        UnsignedInt8* endPtr = dataPtr + numberOfElements * 2;
        while(dataPtr < endPtr) {
          std::swap(*dataPtr, *(dataPtr + 1));
          dataPtr += 2;
        }
        return;
      }


      template <>
      inline void
      genericSwitchByteOrder<4>(UnsignedInt8* dataPtr,
                                size_t numberOfElements,
                                ByteOrder fromByteOrder,
                                ByteOrder toByteOrder)
      {
        if(fromByteOrder == toByteOrder) {
          return;
        }
        UnsignedInt8* endPtr = dataPtr + numberOfElements * 4;
        while(dataPtr < endPtr) {
          std::swap(*dataPtr, *(dataPtr + 3));
          std::swap(*(dataPtr + 1), *(dataPtr + 2));
          dataPtr += 4;
        }
        return;
      }


      template <>
      inline void
      genericSwitchByteOrder<8>(UnsignedInt8* dataPtr,
                                size_t numberOfElements,
                                ByteOrder fromByteOrder,
                                ByteOrder toByteOrder)
      {
        if(fromByteOrder == toByteOrder) {
          return;
        }
        UnsignedInt8* endPtr = dataPtr + numberOfElements * 8;
        while(dataPtr < endPtr) {
          std::swap(*dataPtr, *(dataPtr + 7));
          std::swap(*(dataPtr + 1), *(dataPtr + 6));
          std::swap(*(dataPtr + 2), *(dataPtr + 5));
          std::swap(*(dataPtr + 3), *(dataPtr + 4));
          dataPtr += 8;
        }
        return;
      }


      /**
       * This is an "under-the-hood" function implementation function for
       * the switchByteOrder() function template.  It is not part of the
       * public interface, and may be removed later.
       *
       * The default implementation isn't smart enough to know how to switch
       * any byte ordering.  We use specializations for that.
       **/
      template <int Size>
      inline void
      genericSwitchByteOrder(const UnsignedInt8* fromDataPtr,
                             size_t numberOfElements,
                             UnsignedInt8* toDataPtr,
                             ByteOrder fromByteOrder,
                             ByteOrder toByteOrder)
      {
        if(fromByteOrder == toByteOrder) {
          std::copy(fromDataPtr, fromDataPtr + numberOfElements, toDataPtr);
          return;
        }
        std::ostringstream message;
        message << "This function is not implemented for types of size"
                << Size << ".";
        BRICK_THROW(NotImplementedException, "genericSwitchByteOrder()",
                    message.str().c_str());
      }


      template <>
      inline void
      genericSwitchByteOrder<2>(const UnsignedInt8* fromDataPtr,
                                size_t numberOfElements,
                                UnsignedInt8* toDataPtr,
                                ByteOrder fromByteOrder,
                                ByteOrder toByteOrder)
      {
        if(fromByteOrder == toByteOrder) {
          std::copy(fromDataPtr, fromDataPtr + numberOfElements, toDataPtr);
          return;
        }
        const UnsignedInt8* fromEndPtr =  fromDataPtr + numberOfElements * 2;
        while(fromDataPtr < fromEndPtr) {
          *toDataPtr = *(fromDataPtr + 1);
          *(toDataPtr + 1) = *(fromDataPtr);
          toDataPtr += 2;
          fromDataPtr += 2;
        }
        return;
      }


      template <>
      inline void
      genericSwitchByteOrder<4>(const UnsignedInt8* fromDataPtr,
                                size_t numberOfElements,
                                UnsignedInt8* toDataPtr,
                                ByteOrder fromByteOrder,
                                ByteOrder toByteOrder)
      {
        if(fromByteOrder == toByteOrder) {
          std::copy(fromDataPtr, fromDataPtr + numberOfElements, toDataPtr);
          return;
        }
        const UnsignedInt8* fromEndPtr = fromDataPtr + numberOfElements * 4;
        while(fromDataPtr < fromEndPtr) {
          *toDataPtr = *(fromDataPtr + 3);
          *(toDataPtr + 1) = *(fromDataPtr + 2);
          *(toDataPtr + 2) = *(fromDataPtr + 1);
          *(toDataPtr + 3) = *(fromDataPtr);
          toDataPtr += 4;
          fromDataPtr += 4;
        }
        return;
      }


      template <>
      inline void
      genericSwitchByteOrder<8>(const UnsignedInt8* fromDataPtr,
                                size_t numberOfElements,
                                UnsignedInt8* toDataPtr,
                                ByteOrder fromByteOrder,
                                ByteOrder toByteOrder)
      {
        if(fromByteOrder == toByteOrder) {
          std::copy(fromDataPtr, fromDataPtr + numberOfElements, toDataPtr);
          return;
        }
        const UnsignedInt8* fromEndPtr = fromDataPtr + numberOfElements * 8;
        while(fromDataPtr < fromEndPtr) {
          *toDataPtr = *(fromDataPtr + 7);
          *(toDataPtr + 1) = *(fromDataPtr + 6);
          *(toDataPtr + 2) = *(fromDataPtr + 5);
          *(toDataPtr + 3) = *(fromDataPtr + 4);
          *(toDataPtr + 4) = *(fromDataPtr + 3);
          *(toDataPtr + 5) = *(fromDataPtr + 2);
          *(toDataPtr + 6) = *(fromDataPtr + 1);
          *(toDataPtr + 7) = *(fromDataPtr);
          toDataPtr += 8;
          fromDataPtr += 8;
        }
        return;
      }

    } // namespace privateCode
    /// @endcond

  
    // This function returns the appropriate byte order for the platform
    // on which it is run.
    inline ByteOrder
    getByteOrder()
    {
      UnsignedInt16 byteOrderTester = 0x0102;
      if(*((UnsignedInt8*)(&byteOrderTester)) == 0x01) {
        return BRICK_BIG_ENDIAN;
      }
      return BRICK_LITTLE_ENDIAN;
    }


    // This function swaps the byte order of a single value.
    template <class Type>
    inline Type
    switchByteOrder(Type inputValue,
                    ByteOrder fromByteOrder,
                    ByteOrder toByteOrder)
    {
      Type returnValue;
      switchByteOrder(&inputValue, 1, &returnValue, fromByteOrder, toByteOrder);
      return returnValue;
    }

    
    // This function takes a pointer to a C-style array of values and
    // modifies the array in place so that it has a particular byte
    // order.
    template <class Type>
    void
    switchByteOrder(Type* dataPtr,
                    size_t numberOfElements,
                    ByteOrder fromByteOrder,
                    ByteOrder toByteOrder)
    {
      privateCode::genericSwitchByteOrder<sizeof(Type)>(
        reinterpret_cast<UnsignedInt8*>(dataPtr), numberOfElements,
        fromByteOrder, toByteOrder);
    }

    
    // This function takes a pointer to a C-style array of values and
    // copies it into another C-style array, swapping bytes if necessary
    // so that the output array has the specified byte order.
    template <class Type>
    void
    switchByteOrder(const Type* fromDataPtr,
                    size_t numberOfElements,
                    Type* toDataPtr,
                    ByteOrder fromByteOrder,
                    ByteOrder toByteOrder)
    {
      privateCode::genericSwitchByteOrder<sizeof(Type)>(
        reinterpret_cast<const UnsignedInt8*>(fromDataPtr), numberOfElements,
        reinterpret_cast<UnsignedInt8*>(toDataPtr), fromByteOrder, toByteOrder);
    }

    
  } // namespace common
  
} // namespace brick
    
#endif /* #ifndef BRICK_COMMON_BYTEORDER_HH */
