/**
***************************************************************************
* @file traceableTest.cc
* 
* Source file defining tests for exception trace code.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_USE_TRACEABLE
#define BRICK_COMMON_USE_TRACEABLE 1
#endif /* ifndef BRICK_COMMON_USE_TRACEABLE */

#include <brick/common/exception.hh>
#include <brick/common/traceable.hh>

namespace {

  void
  printDebugInfo(std::string const& traceMessage,
                 std::string const& desiredTraceMessage)
  {
    for(unsigned int ii = 0;
        ii < std::min(traceMessage.size(), desiredTraceMessage.size());
        ++ii) {
      if(traceMessage[ii] != desiredTraceMessage[ii]) {
        std::cout << "Mismatch at " << ii << " "
                  << traceMessage[ii] << " vs. " << desiredTraceMessage[ii]
                  << std::endl;
      }
    }
    if(traceMessage.size() != desiredTraceMessage.size()) {
      std::cout << "Size mismatch: " << traceMessage.size() << ", "
                << desiredTraceMessage.size() << std::endl;
    }
  }

  
  void
  traceFunction(size_t count, const std::string& message, size_t& lineNumber)
  {
    BEGIN_TRACEABLE;
    lineNumber = __LINE__ + 6;
    if(count == 0) {
      BRICK_THROW(brick::common::ValueException, "traceFunction()",
                  message.c_str());
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
    } catch(const brick::common::Exception& caughtException) {
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
        (BRICK_COMMON_TRACE_MESSAGE_LENGTH
         - desiredTraceStream.str().size()
         - ellipsis.size());
      desiredTraceStream << message.substr(0, remainingLength);
      desiredTraceStream << ellipsis;
      // std::cout << caughtException.what() << std::endl;
      // std::string s0 = caughtException.trace();
      // std::string s1 = desiredTraceStream.str();
      // std::cout << s0 << std::endl;
      // std::cout << s1 << std::endl;
      // std::cout << s0.size() << std::endl;
      // std::cout << s1.size() << std::endl;
      std::string traceMessage = brick::common::getTrace(caughtException);
      std::string desiredTraceMessage = desiredTraceStream.str();
      printDebugInfo(traceMessage, desiredTraceMessage);
      return (desiredTraceMessage == traceMessage);
    }
    return false;
  }


  bool
  testLongTrace()
  {
    std::cout << "Testing long trace..." << std::endl;
    
    size_t lineNumber = 0;
    std::ostringstream messageStream;
    for(size_t count = 0; count < BRICK_COMMON_TRACE_MESSAGE_LENGTH;
        ++count) {
      messageStream << (count % 10);
    }
    std::string message = messageStream.str();

    const unsigned int stackLevels =
      BRICK_EXCEPTION_PAYLOAD_SIZE / BRICK_COMMON_TRACE_MESSAGE_LENGTH;
    try {
      traceFunction(stackLevels + 5, message, lineNumber);
    } catch(const brick::common::Exception& caughtException) {
      std::string ellipsis = "...";
      std::string fileName = __FILE__;
      std::ostringstream desiredTraceStream;
      for(size_t count = 0; count < stackLevels; ++count) {
        std::ostringstream subMessageStream;
        subMessageStream
          << "\n"
          << "\n"
          << "  (" << fileName << ", " << lineNumber << "): traceFunction()\n"
          << "    (arg0): " << count << "\n"
          << "    (arg1): ";
        size_t remainingLength = 
          (BRICK_COMMON_TRACE_MESSAGE_LENGTH
           - subMessageStream.str().size()
           - ellipsis.size());
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
      std::string traceMessage = brick::common::getTrace(caughtException);
      std::string desiredTraceMessage = desiredTraceStream.str();
      printDebugInfo(traceMessage, desiredTraceMessage);
      return (desiredTraceMessage == traceMessage);
    }
    return false;
  }


  bool
  testShortTrace()
  {
    std::cout << "Testing short trace..." << std::endl;

    const size_t numLevels = 3;
    size_t lineNumber = 0;
    try {
      traceFunction(numLevels, "-", lineNumber);
    } catch(const brick::common::Exception& caughtException) {
      std::string fileName = __FILE__;
      std::ostringstream desiredTraceStream;
      for(size_t level = 0; level <= numLevels; ++level) {
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
      std::string traceMessage = brick::common::getTrace(caughtException);
      std::string desiredTraceMessage = desiredTraceStream.str();
      printDebugInfo(traceMessage, desiredTraceMessage);
      return (desiredTraceMessage == traceMessage);
    }
    return false;
  }

} // anonymous namespace


// int main(int argc, char** argv)
int main(int, char**)
{
  bool result = true;
  result &= testEllipsis();
  result &= testLongTrace();
  result &= testShortTrace();
  return (result ? 0 : 1);
}
