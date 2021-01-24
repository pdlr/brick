/**
***************************************************************************
* @file brick/common/argumentDescription.hh
*
* Header file declaring functions which can be used to describe
* arguments in stack traces.
*
* Copyright (c) 2005-2007, David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_ARGUMENTDESCRIPTION_HH
#define BRICK_COMMON_ARGUMENTDESCRIPTION_HH

#include <ostream>

namespace brick {

  namespace common {

    /**
     * This function template is responsible deciding how to print
     * argument values during a stack trace.  The default is to use the
     * stream output operator, but you'll want to specialize it for
     * classes that don't implement a stream output operator, and for
     * classes whose stream output operators print inappropriate stuff.
     *
     * @param outputStream This argument is the stream to which the
     * argument description will be sent.
     *
     * @param argument This argument is the actual argument to be described.
     *
     * @return The return value is a reference to parameter outputStream.
     */
    template <class Type>
    std::ostream&
    addArgumentDescription(std::ostream& outputStream, const Type& argument)
    {
      outputStream << argument;
      return outputStream;
    }

  } // namespace common

} // namespace brick


#ifdef BRICK_COMMON_USE_TRACEABLE

#include <map>
#include <vector>

namespace brick {

  namespace common {

    template <class TYPE>
    std::ostream&
    addArgumentDescription(std::ostream& outputStream,
                           const std::vector<TYPE>& argument)
    {
      outputStream << "std::vector([";
      if(argument.empty() == false) {
        outputStream << argument[0];
      }
      for(size_t index0 = 0; index0 < argument.size(); ++index0) {
        outputStream << "," << argument[index0];
      }
      outputStream << "])";
      return outputStream;
    }


    template <class TYPE0, class TYPE1>
    std::ostream&
    addArgumentDescription(std::ostream& outputStream,
                           const std::pair<TYPE0, TYPE1>& argument)
    {
      outputStream << "std::pair(" << argument.first << ", "
                   << argument.second << ")";
      return outputStream;
    }

  } // namespace common

} // namespace brick

#endif /* #ifdef BRICK_COMMON_USE_TRACEABLE */


/* ======= Declarations to maintain compatibility with legacy code. ======= */

namespace brick {

  using common::addArgumentDescription;

} // namespace brick

#endif /* #ifndef BRICK_COMMON_ARGUMENTDESCRIPTION_HH */
