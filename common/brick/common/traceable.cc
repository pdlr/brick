/**
***************************************************************************
* @file brick/common/traceable.cc
*
* Source file defining tools for generating debugging traces when
* failures occur.
*
* Copyright (C) 2005-2011, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/traceable.hh>

namespace brick {

  namespace common {

    // This function is used by BEGIN_TRACEABLE to add
    // levels to an automatically generated stack trace.
    void
    addTrace(brick::common::Exception& caughtException,
             std::string const& message)
    {
      const char* ellipsisString = "...";
      const size_t ellipsisLength = strlen(ellipsisString);

      // Sanity check.
      if(BRICK_COMMON_TRACE_MESSAGE_LENGTH < ellipsisLength) {
        return;
      }

      // Figure out if we need to add ellipsis.
      bool useEllipsis = false;
      size_t messageLength = message.size();
      if(messageLength > BRICK_COMMON_TRACE_MESSAGE_LENGTH) {
        useEllipsis = true;
        messageLength = BRICK_COMMON_TRACE_MESSAGE_LENGTH - ellipsisLength;
      }

      // Copy the message into the next part of the stack trace.
      // Exception is smart enough not to overflow.
      caughtException.setPayload(
        caughtException.getPayloadSize(), message.c_str(), messageLength);
      if(useEllipsis) {
        caughtException.setPayload(
          caughtException.getPayloadSize(), ellipsisString, ellipsisLength);
      }
    }


    // This function returns the stack trace accumulated by
    // BEGIN_TRACEABLE/END_TRACEABLE.
    std::string
    getTrace(brick::common::Exception const& caughtException)
    {
      unsigned int dummy;
      char* buffer = new char[caughtException.getPayloadSize()];
      try {
        caughtException.getPayload(buffer, dummy);
        std::string returnValue(buffer, caughtException.getPayloadSize());
        delete[] buffer;
        return returnValue;
      } catch(...) {
        delete[] buffer;
        throw;
      }
      return "";
    }

  } // namespace common

} // namespace brick
