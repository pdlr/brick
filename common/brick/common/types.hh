/**
***************************************************************************
* @file brick/common/types.hh
*
* Header file defining explicit types, such as Int32, Int64, etc.
*
* Copyright (C) 2005, 2010, David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_TYPES_HH
#define BRICK_COMMON_TYPES_HH

namespace brick {

  namespace common {

    /// @cond privateCode
    namespace privateCode {

      // This struct template must be specialized once for each
      // architecture under which we need to build in order to define
      // the sizes of numeric types.  Except for these
      // specializations, GlobalTypeStruct is never used directly.
      template<int TMPL_SIZEOF_CHAR,
               int TMPL_SIZEOF_SHORT_INT,
               int TMPL_SIZEOF_INT,
               int TMPL_SIZEOF_LONG_INT,
               int TMPL_SIZEOF_LONG_LONG_INT,
               int TMPL_SIZEOF_FLOAT,
               int TMPL_SIZEOF_DOUBLE>
      struct GlobalTypeStruct {};


      // This template specialization defines type names for 32 bit
      // i386 machines.
      template<>
      struct GlobalTypeStruct<1, 2, 4, 4, 8, 4, 8> {
        typedef char                    Int8;
        typedef short int               Int16;
        typedef int                     Int32;
        typedef long long int           Int64;
        typedef float                   Float32;
        typedef double                  Float64;
        typedef unsigned char           UnsignedInt8;
        typedef unsigned short int      UnsignedInt16;
        typedef unsigned int            UnsignedInt32;
        typedef unsigned long long int  UnsignedInt64;
      };


      // This template specialization defines type names for x86_64
      // machines, such as AMD64 and 64-bit Intel Xeon.
      template<>
      struct GlobalTypeStruct<1, 2, 4, 8, 8, 4, 8> {
        typedef char                    Int8;
        typedef short int               Int16;
        typedef int                     Int32;
        typedef long int                Int64;
        typedef float                   Float32;
        typedef double                  Float64;
        typedef unsigned char           UnsignedInt8;
        typedef unsigned short int      UnsignedInt16;
        typedef unsigned int            UnsignedInt32;
        typedef unsigned long int       UnsignedInt64;
      };


      // This typedef just makes it easier to write the whole slew of
      // upcoming typedefs.
      typedef
      GlobalTypeStruct<sizeof(char),
                       sizeof(short int),
                       sizeof(int),
                       sizeof(long int),
                       sizeof(long long int),
                       sizeof(float),
                       sizeof(double)>
      LocalTypeStruct;

    } // namespace privateCode
    /// @endcond


    // These are the typedefs we'll actually use in the rest of our
    // code.
    typedef privateCode::LocalTypeStruct::Int8 Int8;
    typedef privateCode::LocalTypeStruct::Int16 Int16;
    typedef privateCode::LocalTypeStruct::Int32 Int32;
    typedef privateCode::LocalTypeStruct::Int64 Int64;
    typedef privateCode::LocalTypeStruct::Float32 Float32;
    typedef privateCode::LocalTypeStruct::Float64 Float64;
    typedef privateCode::LocalTypeStruct::UnsignedInt8 UnsignedInt8;
    typedef privateCode::LocalTypeStruct::UnsignedInt16 UnsignedInt16;
    typedef privateCode::LocalTypeStruct::UnsignedInt32 UnsignedInt32;
    typedef privateCode::LocalTypeStruct::UnsignedInt64 UnsignedInt64;

    // Shorter type names for more convenient use.
    typedef privateCode::LocalTypeStruct::UnsignedInt8 UInt8;
    typedef privateCode::LocalTypeStruct::UnsignedInt16 UInt16;
    typedef privateCode::LocalTypeStruct::UnsignedInt32 UInt32;
    typedef privateCode::LocalTypeStruct::UnsignedInt64 UInt64;

  } // namespace common

}  // namespace brick

#endif /* #ifndef BRICK_COMMON_TYPES_HH */
