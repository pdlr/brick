/**
***************************************************************************
* @file traceableTest.cc
* Source file defining tests for exception trace code.
*
* Copyright (C) 2005 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_USE_TRACEABLE
#define BRICK_COMMON_USE_TRACEABLE
#endif /* ifndef BRICK_COMMON_USE_TRACEABLE */

#include <brick/common/exception.hh>
#include <brick/common/traceable.hh>

namespace brick {

  void
  traceFunction(size_t count, const std::string& message, size_t& lineNumber)
  {
    BEGIN_TRACEABLE;
    lineNumber = __LINE__ + 5;
    if(count == 0) {
      BRICK_THROW(ValueException, "traceFunction()", message.c_str());
    }
    traceFunction(count - 1, message, lineNumber);
    END_TRACEABLE(traceFunction, (count, message, lineNumber));
  }


  bool
  testEllipsis()
  {
    std::cout << "Testing ellipsis..." << std::endl;
    
    size_t lineNumber = 0;
    std::ostringstream messageStream;
    for(size_t count = 0; count < BRICK_COMMON_TRACE_MESSAGE_LENGTH;
        ++count) {
      messageStream << (count % 10);
    }
    std::string message = messageStream.str();

    try {
      traceFunction(0, message, lineNumber);
    } catch(const brick::Exception& caughtException) {
      std::string ellipsis = "...";
      std::string fileName = __FILE__;
      std::ostringstream desiredTraceStream;
      desiredTraceStream
        << "\n"
        << "\n"
        << "  (" << fileName << ", " << lineNumber << "): traceFunction()\n"
        << "    (arg0): " << 0 << "\n"
        << "    (arg1): ";
      size_t remainingLength =
        (BRICK_EXCEPTION_TRACE_MESSAGE_LENGTH
         - desiredTraceStream.str().size()
         - ellipsis.size()
         - 1);
      desiredTraceStream << message.substr(0, remainingLength);
      desiredTraceStream << ellipsis;
      // std::cout << caughtException.what() << std::endl;
      // std::string s0 = caughtException.trace();
      // std::string s1 = desiredTraceStream.str();
      // std::cout << s0 << std::endl;
      // std::cout << s1 << std::endl;
      // std::cout << s0.size() << std::endl;
      // std::cout << s1.size() << std::endl;
      return desiredTraceStream.str() == caughtException.trace();
    }
    return false;
  }


  bool
  testLongTrace()
  {
    std::cout << "Testing long trace..." << std::endl;
    
    size_t lineNumber = 0;
    std::ostringstream messageStream;
    for(size_t count = 0; count < BRICK_EXCEPTION_TRACE_MESSAGE_LENGTH;
        ++count) {
      messageStream << (count % 10);
    }
    std::string message = messageStream.str();

    try {
      traceFunction(BRICK_EXCEPTION_TRACE_REQUIRED_STACK_LEVELS + 5,
                    message, lineNumber);
    } catch(const brick::Exception& caughtException) {
      std::string ellipsis = "...";
      std::string fileName = __FILE__;
      std::ostringstream desiredTraceStream;
      for(size_t count = 0;
          count < BRICK_EXCEPTION_TRACE_REQUIRED_STACK_LEVELS;
          ++count) {
        std::ostringstream subMessageStream;
        subMessageStream
          << "\n"
          << "\n"
          << "  (" << fileName << ", " << lineNumber << "): traceFunction()\n"
          << "    (arg0): " << count << "\n"
          << "    (arg1): ";
        size_t remainingLength = 
          (BRICK_EXCEPTION_TRACE_MESSAGE_LENGTH
           - subMessageStream.str().size()
           - ellipsis.size()
           - 1);
        subMessageStream << message.substr(0, remainingLength);
        subMessageStream << ellipsis;
        desiredTraceStream << subMessageStream.str();
      }
      // std::cout << caughtException.what() << std::endl;
      // std::string s0 = caughtException.trace();
      // std::string s1 = desiredTraceStream.str();
      // std::cout << s0 << std::endl;
      // std::cout << s1 << std::endl;
      // std::cout << s0.size() << std::endl;
      // std::cout << s1.size() << std::endl;
      return desiredTraceStream.str() == caughtException.trace();
    }
    return false;
  }


  bool
  testShortTrace()
  {
    std::cout << "Testing short trace..." << std::endl;
    
    size_t lineNumber = 0;
    try {
      traceFunction(3, "-", lineNumber);
    } catch(const brick::Exception& caughtException) {
      std::string fileName = __FILE__;
      std::ostringstream desiredTraceStream;
      for(size_t level = 0; level <= 3; ++level) {
        desiredTraceStream
          << "\n"
          << "\n"
          << "  (" << fileName << ", " << lineNumber << "): traceFunction()\n"
          << "    (arg0): " << level << "\n"
          << "    (arg1): -\n"
          << "    (arg2): " << lineNumber;
      }
      // std::cout << caughtException.what() << std::endl;
      // std::cout << caughtException.trace() << std::endl;
      // std::cout << desiredTraceStream.str() << std::endl;
      return desiredTraceStream.str() == caughtException.trace();
    }
    return false;
  }

} // namespace brick


int main(int argc, char** argv)
{
  bool result = true;
  result &= brick::testEllipsis();
  result &= brick::testLongTrace();
  result &= brick::testShortTrace();
  return (result ? 0 : 1);
}
