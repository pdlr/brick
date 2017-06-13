/**
***************************************************************************
* @file brick/numeric/numericTraits.hh
*
* Header file declaring traits classes for handling type promotion and
* similar tasks.
*
* Copyright (C) 2003-2017, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_NUMERICTRAITS_HH
#define BRICK_NUMERIC_NUMERICTRAITS_HH

#include <limits>

#include <brick/common/exception.hh>
#include <brick/common/types.hh>

namespace brick {

  namespace numeric {

    /**
     ** This class template allows generic code to be written that
     ** handles arithmetic operations without overflowing or losing
     ** precision.  It provides typedefs for controlling the precision
     ** of computations, numeric limits, etc.  It is specialized
     ** explicitly for several numeric types below.
     **/
    template <class Type0, class Type1>
    struct ArithmeticTraits {

      /**
       * This typedef is used to control the precision of the result
       * when values of Type0 and Type1 are multiplied.  ProductType
       * should have enough precision that the result of multiplying a
       * Type0 instance by a Type1 instance can be represented without
       * overflow or loss of information.
       */
      typedef Type0 ProductType;

      
      /**
       * This typedef is used to control the precision of the result
       * when values of Type0 and Type1 are added.  SumType should
       * have enough precision that the result of adding a Type0
       * instance to a Type1 instance can be represented without
       * overflow or loss of information.
       */
      typedef Type0 SumType;
    };


    /**
     ** This is a base class for NumericTraits, and should generally
     ** not be used directly in client code.
     **/
    template <class Type>
    struct NumericTraitsBase {

      /**
       * This typedef is used to control the way arrays of Type are
       * printed by operator<<().  For example, when a char is
       * serialized for stream output as part of an array, each char
       * is formatted as a number, rather than as an ascii symbol.
       * Illustration -- you might see:
       *
       * @code
       * Array1D([64, 97, 32])
       * @endcode
       *
       * rather than:
       *
       * @code
       *   Array1D([@, a,  ])
       * @endcode
       */
      typedef Type TextOutputType;


      /** 
       * Returns the a very small number, the positive difference
       * between 1 and the next biggest representable value.  Use this
       * in place of FLT_EPSILON, DBL_EPSILON, etc.
       * 
       * @return The return value is the epsilon value.
       */
      static inline Type
      epsilon() {
        return std::numeric_limits<Type>::epsilon();
      }
      

      /** 
       * This member function indicates whether the specified type is an
       * integer type or not.  It is not implemented for the general case,
       * so NumericTraits must be specialized for every class which requires
       * this functionality.
       * 
       * @return The return value should true if this type is an integer
       * type, false other wise.  The general version of this function just
       * throws and exception.
       */
      static inline bool
      isIntegral() {
        BRICK_THROW(common::NotImplementedException,
                    "NumericTraits::isIntegral()",
                    "In order to use isIntegral, NumericTraits must be "
                    "specialized for each type.");
        return true; // Appease the compiler.
      }

    };
    
    /**
     ** This class template allows generic numerical code to be
     ** written by providing information about the numeric type of its
     ** template argument.  It is specialized explicitly for several
     ** numeric types below.
     **/
    template <class Type>
    struct NumericTraits : public NumericTraitsBase<Type> {
      // Specializations will override inherited members from
      // NumericTraitsBase.  See specializations below for examples.
    };


    // == Specializations of ArithmeticTraits == 
    template <>
    struct ArithmeticTraits<common::Int8, common::Int8> {
      typedef common::Int16 ProductType;
      typedef common::Int16 SumType;
    };


    template <>
    struct ArithmeticTraits<common::UnsignedInt8, common::UnsignedInt8> {
      typedef common::UnsignedInt16 ProductType;
      typedef common::UnsignedInt16 SumType;
    };


    template <>
    struct ArithmeticTraits<common::Int16, common::Int16> {
      typedef common::Int32 ProductType;
      typedef common::Int32 SumType;
    };


    template <>
    struct ArithmeticTraits<common::UnsignedInt16, common::UnsignedInt16> {
      typedef common::UnsignedInt32 ProductType;
      typedef common::UnsignedInt32 SumType;
    };


    template <>
    struct ArithmeticTraits<common::Int32, common::Int32> {
      typedef common::Int64 ProductType;
      typedef common::Int64 SumType;
    };
    
    
    template <>
    struct ArithmeticTraits<common::UnsignedInt32, common::UnsignedInt32> {
      typedef common::UnsignedInt64 ProductType;
      typedef common::UnsignedInt64 SumType;
    };
    
    
    template <>
    struct ArithmeticTraits<common::Int64, common::Int64> {
      typedef common::Int64 ProductType;
      typedef common::Int64 SumType;
    };
    
    
    template <>
    struct ArithmeticTraits<common::UnsignedInt64, common::UnsignedInt64> {
      typedef common::UnsignedInt64 ProductType;
      typedef common::UnsignedInt64 SumType;
    };
    
    
    template <>
    struct ArithmeticTraits<common::Float32, common::Float32> {
      typedef common::Float32 ProductType;
      typedef common::Float32 SumType;
    };
    
    
    template <>
    struct ArithmeticTraits<common::Float64, common::Float64> {
      typedef common::Float64 ProductType;
      typedef common::Float64 SumType;
    };
    
    
    // == Specializations of NumericTriats == 

    template <>
    struct NumericTraits<char>
      : public NumericTraitsBase<char>
    {
      typedef int TextOutputType;
      static inline bool isIntegral() {return true;}
    };

    
    template <>
    struct NumericTraits<unsigned char>
      : public NumericTraitsBase<unsigned char>
    {
      typedef unsigned int TextOutputType;
      static inline bool isIntegral() {return true;}
    };

  
    template <>
    struct NumericTraits<common::Int16>
      : public NumericTraitsBase<common::Int16>
    {
      static inline bool isIntegral() {return true;}
    };


    template <>
    struct NumericTraits<common::UnsignedInt16>
      : public NumericTraitsBase<common::UnsignedInt16>
    {
      static inline bool isIntegral() {return true;}
    };


    template <>
    struct NumericTraits<common::Int32>
      : public NumericTraitsBase<common::Int32>
    {
      static inline bool isIntegral() {return true;}
    };


    template <>
    struct NumericTraits<common::UnsignedInt32>
      : public NumericTraitsBase<common::UnsignedInt32>
    {
      static inline bool isIntegral() {return true;}
    };


    template <>
    struct NumericTraits<common::Int64>
      : public NumericTraitsBase<common::Int64>
    {
      static inline bool isIntegral() {return true;}
    };


    template <>
    struct NumericTraits<common::UnsignedInt64>
      : public NumericTraitsBase<common::UnsignedInt64>
    {
      static inline bool isIntegral() {return true;}
    };


    template <>
    struct NumericTraits<common::Float32>
      : public NumericTraitsBase<common::Float32>
    {
      static inline bool isIntegral() {return false;}
    };


    template <>
    struct NumericTraits<common::Float64>
      : public NumericTraitsBase<common::Float64>
    {
      static inline bool isIntegral() {return false;}
    };

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_NUMERICTRAITS_HH */
