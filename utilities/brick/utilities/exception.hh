/**
***************************************************************************
* @file brick/utilities/exception.hh
*
* Header file declaring some exception types and providing namespace
* documentation.
*
* Copyright (c) 2007, David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_EXCEPTION_HH
#define BRICK_UTILITIES_EXCEPTION_HH

#include <brick/common/exception.hh>

namespace brick {

  /**
   ** This namespace contains general utility functions and classes
   ** for string manipluation, program option parsing, working with
   ** time and filenames, and a few other odds and ends.
   **/
  namespace utilities {

    /**
     ** This is an Exception class for errors involving failed
     ** conversions from one type to another.
     **/
    class ConversionException;
    BRICK_DECLARE_EXCEPTION_TYPE(
      ConversionException, brick::common::ValueException);

  } // namespace utilities

} // namespace brick

#endif /* #ifndef BRICK_UTILITIES_EXCEPTION_HH */
