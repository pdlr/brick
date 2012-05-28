/**
***************************************************************************
* @file brick/numeric/typeConverter.hh
*
* Header file declaring TypeConverter class template.
*
* Copyright (C) 2012, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_TYPECONVERTER_HH
#define BRICK_NUMERIC_TYPECONVERTER_HH

#include <brick/common/types.hh>

namespace brick {

  namespace numeric {

    /**
     ** This class template allows generic code to be written that
     ** handles conversions between numeric types.  It is specialized
     ** explicitly for several numeric types below.
     **/
    template <class Type0, class Type1>
    struct TypeConverter {
      Type0 operator()(Type1 const& arg) {return static_cast<Type0>(arg);}
    };


    /** 
     * Convenience function to allow syntactically pleasant use of the
     * TypeConverter class template.
     * 
     * @param arg This argument is the value to be converted.
     * 
     * @return The return value is the result of the conversion.
     */
    template <class Type0, class Type1>
    Type0
    convertType(Type1 const& arg) {
      return (TypeConverter<Type0, Type1>()).operator()(arg);
    }

    
    // == Specializations of TypeConverter == 

    template <>
    struct TypeConverter<common::Int8, common::Float32> {
      common::Int8 operator()(common::Float32 const& arg) {
        return static_cast<common::Int8>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::Int8, common::Float64> {
      common::Int8 operator()(common::Float32 const& arg) {
        return static_cast<common::Int8>(arg + 0.5);
      }
    };
    
    template <>
    struct TypeConverter<common::Int16, common::Float32> {
      common::Int8 operator()(common::Float32 const& arg) {
        return static_cast<common::Int16>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::Int16, common::Float64> {
      common::Int8 operator()(common::Float32 const& arg) {
        return static_cast<common::Int16>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::Int32, common::Float32> {
      common::Int8 operator()(common::Float32 const& arg) {
        return static_cast<common::Int32>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::Int32, common::Float64> {
      common::Int8 operator()(common::Float32 const& arg) {
        return static_cast<common::Int32>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::Int64, common::Float32> {
      common::Int8 operator()(common::Float32 const& arg) {
        return static_cast<common::Int64>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::Int64, common::Float64> {
      common::Int8 operator()(common::Float32 const& arg) {
        return static_cast<common::Int64>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::UInt8, common::Float32> {
      common::UInt8 operator()(common::Float32 const& arg) {
        return static_cast<common::UInt8>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::UInt8, common::Float64> {
      common::UInt8 operator()(common::Float32 const& arg) {
        return static_cast<common::UInt8>(arg + 0.5);
      }
    };
    
    template <>
    struct TypeConverter<common::UInt16, common::Float32> {
      common::UInt8 operator()(common::Float32 const& arg) {
        return static_cast<common::UInt16>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::UInt16, common::Float64> {
      common::UInt8 operator()(common::Float32 const& arg) {
        return static_cast<common::UInt16>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::UInt32, common::Float32> {
      common::UInt8 operator()(common::Float32 const& arg) {
        return static_cast<common::UInt32>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::UInt32, common::Float64> {
      common::UInt8 operator()(common::Float32 const& arg) {
        return static_cast<common::UInt32>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::UInt64, common::Float32> {
      common::UInt8 operator()(common::Float32 const& arg) {
        return static_cast<common::UInt64>(arg + 0.5);
      }
    };

    template <>
    struct TypeConverter<common::UInt64, common::Float64> {
      common::UInt8 operator()(common::Float32 const& arg) {
        return static_cast<common::UInt64>(arg + 0.5);
      }
    };

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_TYPECONVERTER_HH */
