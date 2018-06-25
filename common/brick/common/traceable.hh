/**
***************************************************************************
* @file brick/common/traceable.hh
*
* Header file declaring tools for generating debugging traces when
* failures occur.
*
* Copyright (C) 2005-2011, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_TRACEABLE_HH
#define BRICK_COMMON_TRACEABLE_HH

#include <iostream>
#include <sstream>
#include <vector>
#include <brick/common/argumentDescription.hh>
#include <brick/common/exception.hh>

// This file defines some preprocessor macros.  The definitions of
// these macros depend on whether or not we want to generate traces.
// We default to zero.


#ifndef BRICK_COMMON_USE_TRACEABLE
/**
 ** This macro controls whether or not to the BEGIN_TRACEABLE and
 ** END_TRACEABLE macros should be implemented.
 **
 ** WARNING: This functionality uses the setPayload()/getPayload()
 ** feature of brick::common::Exception.  If you are already using the
 ** payload in some other way, you should leave
 ** BRICK_COMMON_USE_TRACEABLE at its default value of 0.
 **/
#define BRICK_COMMON_USE_TRACEABLE 0
#endif /* #ifndef BRICK_COMMON_USE_TRACEABLE */


/**
 ** This macro specifies how many characters are allowable in the
 ** description of each level of a stack trace.  The total number of
 ** stack levels you can report will be greater than or equal to
 ** BRICK_EXCEPTION_PAYLOAD_SIZE / BRICK_COMMON_TRACE_MESSAGE_LENGTH.
 **/
#define BRICK_COMMON_TRACE_MESSAGE_LENGTH 512


/**
 ** If BRICK_COMMON_USE_TRACEABLE is nonzero, this macro works with
 ** END_TRACEABLE to build a stack trace.  Use it inside each function
 ** definition, like this:
 **
 **   int myFunction(int argument0, int argument1)
 **   {
 **     BEGIN_TRACEABLE;
 **     // Your code goes here.
 **     return 0;
 **     END_TRACEABLE(myFunction, (argument0, argument1));
 **   }
 **/
#if BRICK_COMMON_USE_TRACEABLE

#define BEGIN_TRACEABLE \
  try {

#else /* #if BRICK_COMMON_USE_TRACEABLE */

#define BEGIN_TRACEABLE

#endif /* #if BRICK_COMMON_USE_TRACEABLE ... #else */

/**
 ** If BRICK_COMMON_USE_TRACEABLE is defined, this macro works with
 ** BEGIN_TRACEABLE to build the stack trace.  Use it inside each
 ** function definition, like this:
 **
 ** @code
 **   int myFunction(int argument0, int argument1)
 **   {
 **     BEGIN_TRACEABLE;
 **     // Your code goes here.
 **     return 0;
 **     END_TRACEABLE(myFunction, (argument0, argument1));
 **   }
 **
 **   try {
 **     myFunction(argument0, argument1);
 **   } catch(brick::common::Exception& caughtException) {
 **     // Print the trace message.
 **     std::cerr << brick::common::getTrace(caughtException) << std::endl;
 **   }
 ** @endcode
 **
 ** WARNING: This functionality uses the setPayload()/getPayload()
 ** feature of brick::common::Exception.  If you are already using the
 ** payload in some other way, you cannot use the BEGIN_TRACEABLE.
 **/
#if BRICK_COMMON_USE_TRACEABLE

#define END_TRACEABLE(functionName, argumentList)                          \
  } catch(brick::common::Exception& caughtException) {                     \
    try {                                                                  \
      std::ostringstream brick_common_end_traceable_message;               \
      brick_common_end_traceable_message                                   \
        << "\n\n  (" << __FILE__ << ", " << __LINE__ << "): "              \
        << #functionName << "()";                                          \
      brick_common_end_traceable_message                                   \
        << brick::common::describeArguments argumentList;                  \
      addTrace(caughtException, brick_common_end_traceable_message.str()); \
    } catch(...) {                                                         \
      /* Empty. The most likely reason to get here is std::bad_alloc   */  \
      /* during the call to describeArguments() or during the call to  */  \
      /* addTrace().                                                   */  \
    }                                                                      \
    /* This throw statment should rethrow the brick::common Exception  */  \
    /* instance, _not_ any exception caught by the catch(...)          */  \
    /* statement above.                                                */  \
    throw;                                                                 \
  }

#else /* #if BRICK_COMMON_USE_TRACEABLE */

#define END_TRACEABLE(functionName, argumentList)

#endif /* #if BRICK_COMMON_USE_TRACEABLE ... #else */

namespace brick {

  namespace common {

    /**
     * This function is used by BEGIN_TRACEABLE to add
     * levels to an automatically generated stack trace.
     *
     * WARNING: This function uses the setPayload()/getPayload()
     * feature of brick::common::Exception.  If you are already using
     * the payload in some other way, you cannot use addTrace().
     *
     * @param caughtException This argument is the exception whose
     * payload is to be updated.
     *
     * @param traceMessage This argument is the text to be added to
     * the stack trace.
     */
    void
    addTrace(brick::common::Exception& caughtException,
             std::string const& message);


    /**
     * This function returns the stack trace accumulated by
     * BEGIN_TRACEABLE/END_TRACEABLE.
     *
     * @param caughtException This argument is the exception from
     * which to read the trace.
     *
     * @return The return value is a string representation of the
     * stack trace.
     */
    std::string
    getTrace(brick::common::Exception const& caughtException);


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @return Since this function has no arguments, the return value is
     * an empty string.
     */
    inline std::string
    describeArguments()
    {
      return std::string("");
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0>
    std::string
    describeArguments(const Type0& argument0)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument11 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10, class Type11>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10,
		      const Type11& argument11)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      message << "\n    (arg11): ";
      addArgumentDescription(message, argument11);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument11 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument12 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10, class Type11, class Type12>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10,
		      const Type11& argument11,
		      const Type12& argument12)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      message << "\n    (arg11): ";
      addArgumentDescription(message, argument11);
      message << "\n    (arg12): ";
      addArgumentDescription(message, argument12);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument11 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument12 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument13 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10, class Type11, class Type12, class Type13>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10,
		      const Type11& argument11,
		      const Type12& argument12,
		      const Type13& argument13)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      message << "\n    (arg11): ";
      addArgumentDescription(message, argument11);
      message << "\n    (arg12): ";
      addArgumentDescription(message, argument12);
      message << "\n    (arg13): ";
      addArgumentDescription(message, argument13);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument11 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument12 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument13 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument14 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10, class Type11, class Type12, class Type13, class Type14>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10,
		      const Type11& argument11,
		      const Type12& argument12,
		      const Type13& argument13,
		      const Type14& argument14)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      message << "\n    (arg11): ";
      addArgumentDescription(message, argument11);
      message << "\n    (arg12): ";
      addArgumentDescription(message, argument12);
      message << "\n    (arg13): ";
      addArgumentDescription(message, argument13);
      message << "\n    (arg14): ";
      addArgumentDescription(message, argument14);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument11 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument12 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument13 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument14 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument15 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10, class Type11, class Type12, class Type13, class Type14,
              class Type15>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10,
		      const Type11& argument11,
		      const Type12& argument12,
		      const Type13& argument13,
		      const Type14& argument14,
		      const Type15& argument15)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      message << "\n    (arg11): ";
      addArgumentDescription(message, argument11);
      message << "\n    (arg12): ";
      addArgumentDescription(message, argument12);
      message << "\n    (arg13): ";
      addArgumentDescription(message, argument13);
      message << "\n    (arg14): ";
      addArgumentDescription(message, argument14);
      message << "\n    (arg15): ";
      addArgumentDescription(message, argument15);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument11 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument12 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument13 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument14 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument15 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument16 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10, class Type11, class Type12, class Type13, class Type14,
              class Type15, class Type16>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10,
		      const Type11& argument11,
		      const Type12& argument12,
		      const Type13& argument13,
		      const Type14& argument14,
		      const Type15& argument15,
		      const Type16& argument16)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      message << "\n    (arg11): ";
      addArgumentDescription(message, argument11);
      message << "\n    (arg12): ";
      addArgumentDescription(message, argument12);
      message << "\n    (arg13): ";
      addArgumentDescription(message, argument13);
      message << "\n    (arg14): ";
      addArgumentDescription(message, argument14);
      message << "\n    (arg15): ";
      addArgumentDescription(message, argument15);
      message << "\n    (arg16): ";
      addArgumentDescription(message, argument16);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument11 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument12 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument13 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument14 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument15 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument16 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument17 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10, class Type11, class Type12, class Type13, class Type14,
              class Type15, class Type16, class Type17>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10,
		      const Type11& argument11,
		      const Type12& argument12,
		      const Type13& argument13,
		      const Type14& argument14,
		      const Type15& argument15,
		      const Type16& argument16,
		      const Type17& argument17)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      message << "\n    (arg11): ";
      addArgumentDescription(message, argument11);
      message << "\n    (arg12): ";
      addArgumentDescription(message, argument12);
      message << "\n    (arg13): ";
      addArgumentDescription(message, argument13);
      message << "\n    (arg14): ";
      addArgumentDescription(message, argument14);
      message << "\n    (arg15): ";
      addArgumentDescription(message, argument15);
      message << "\n    (arg16): ";
      addArgumentDescription(message, argument16);
      message << "\n    (arg17): ";
      addArgumentDescription(message, argument17);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument11 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument12 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument13 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument14 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument15 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument16 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument17 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument18 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10, class Type11, class Type12, class Type13, class Type14,
              class Type15, class Type16, class Type17, class Type18>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10,
		      const Type11& argument11,
		      const Type12& argument12,
		      const Type13& argument13,
		      const Type14& argument14,
		      const Type15& argument15,
		      const Type16& argument16,
		      const Type17& argument17,
		      const Type18& argument18)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      message << "\n    (arg11): ";
      addArgumentDescription(message, argument11);
      message << "\n    (arg12): ";
      addArgumentDescription(message, argument12);
      message << "\n    (arg13): ";
      addArgumentDescription(message, argument13);
      message << "\n    (arg14): ";
      addArgumentDescription(message, argument14);
      message << "\n    (arg15): ";
      addArgumentDescription(message, argument15);
      message << "\n    (arg16): ";
      addArgumentDescription(message, argument16);
      message << "\n    (arg17): ";
      addArgumentDescription(message, argument17);
      message << "\n    (arg18): ";
      addArgumentDescription(message, argument18);
      return message.str();
    }


    /**
     * This function template formats its entire argument list for
     * inclusion in a stack trace.  You'll need to have templates like
     * this defined which take 1, 2, 3, 4, etc., arguments.
     *
     * @param argument0 This argument is the first argument to be included
     * in the formatted output.
     *
     * @param argument1 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument2 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument3 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument4 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument5 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument6 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument7 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument8 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument9 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument10 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument11 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument12 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument13 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument14 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument15 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument16 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument17 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument18 This argument is the second argument to be included
     * in the formatted output.
     *
     * @param argument19 This argument is the second argument to be included
     * in the formatted output.
     *
     * @return The return value is a string describing the argument list.
     */
    template <class Type0, class Type1, class Type2, class Type3, class Type4,
              class Type5, class Type6, class Type7, class Type8, class Type9,
              class Type10, class Type11, class Type12, class Type13, class Type14,
              class Type15, class Type16, class Type17, class Type18, class Type19>
    std::string
    describeArguments(const Type0& argument0,
		      const Type1& argument1,
		      const Type2& argument2,
		      const Type3& argument3,
		      const Type4& argument4,
		      const Type5& argument5,
		      const Type6& argument6,
		      const Type7& argument7,
		      const Type8& argument8,
		      const Type9& argument9,
		      const Type10& argument10,
		      const Type11& argument11,
		      const Type12& argument12,
		      const Type13& argument13,
		      const Type14& argument14,
		      const Type15& argument15,
		      const Type16& argument16,
		      const Type17& argument17,
		      const Type18& argument18,
		      const Type19& argument19)
    {
      std::ostringstream message;
      message << "\n    (arg0): ";
      addArgumentDescription(message, argument0);
      message << "\n    (arg1): ";
      addArgumentDescription(message, argument1);
      message << "\n    (arg2): ";
      addArgumentDescription(message, argument2);
      message << "\n    (arg3): ";
      addArgumentDescription(message, argument3);
      message << "\n    (arg4): ";
      addArgumentDescription(message, argument4);
      message << "\n    (arg5): ";
      addArgumentDescription(message, argument5);
      message << "\n    (arg6): ";
      addArgumentDescription(message, argument6);
      message << "\n    (arg7): ";
      addArgumentDescription(message, argument7);
      message << "\n    (arg8): ";
      addArgumentDescription(message, argument8);
      message << "\n    (arg9): ";
      addArgumentDescription(message, argument9);
      message << "\n    (arg10): ";
      addArgumentDescription(message, argument10);
      message << "\n    (arg11): ";
      addArgumentDescription(message, argument11);
      message << "\n    (arg12): ";
      addArgumentDescription(message, argument12);
      message << "\n    (arg13): ";
      addArgumentDescription(message, argument13);
      message << "\n    (arg14): ";
      addArgumentDescription(message, argument14);
      message << "\n    (arg15): ";
      addArgumentDescription(message, argument15);
      message << "\n    (arg16): ";
      addArgumentDescription(message, argument16);
      message << "\n    (arg17): ";
      addArgumentDescription(message, argument17);
      message << "\n    (arg18): ";
      addArgumentDescription(message, argument18);
      message << "\n    (arg19): ";
      addArgumentDescription(message, argument19);
      return message.str();
    }

  } // namespace common

} // namespace brick

#endif /* #ifndef BRICK_COMMON_TRACEABLE_HH */
