/**
***************************************************************************
* @file brick/common/exception.cc
* 
* Source file defining some exception types.
*
* Copyright (C) 2003-2011, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <brick/common/exception.hh>

// Anonymous namespace for local functions.
namespace {

  // std::snprintf() doesn't seem to be supported on some
  // architectures, so we define a crippled version here which fills
  // our needs.
  int l_snprintf(char* targetPtr, size_t targetSize, 
                 const char* formatPtr, ...) {
    char workingBuffer[1024];
    va_list varArgsList;
    va_start(varArgsList, formatPtr);

    const char* workingFormatPtr = formatPtr;
    char* workingTargetPtr = targetPtr;
    size_t spaceLeft = targetSize - 1;

    // Stop if we've used up either the available space or the format string.
    while(spaceLeft && (*workingFormatPtr != '\0')) {
      int intValue;
      double doubleValue;
      char* charPtrValue;
      size_t stringSize;
      if(*workingFormatPtr == '%') {
        ++workingFormatPtr;
        switch(*workingFormatPtr) {
        case '\0':
          break;
        case 'd':
          intValue = va_arg(varArgsList, int);
          sprintf(workingBuffer, "%d", intValue);
          std::strncpy(workingTargetPtr, workingBuffer, spaceLeft);
          stringSize = strlen(workingBuffer);
          workingTargetPtr +=
            (stringSize > spaceLeft) ? spaceLeft : stringSize;
          spaceLeft -= (stringSize > spaceLeft) ? spaceLeft : stringSize;
          ++workingFormatPtr;
          break;
        case 'f':
          doubleValue = va_arg(varArgsList, double);
          sprintf(workingBuffer, "%f", doubleValue);
          std::strncpy(workingTargetPtr, workingBuffer, spaceLeft);
          stringSize = strlen(workingBuffer);
          workingTargetPtr +=
            (stringSize > spaceLeft) ? spaceLeft : stringSize;
          spaceLeft -= (stringSize > spaceLeft) ? spaceLeft : stringSize;
          ++workingFormatPtr;
          break;
        case 's':
          charPtrValue = va_arg(varArgsList, char*);
          std::strncpy(workingTargetPtr, charPtrValue, spaceLeft);
          stringSize = strlen(charPtrValue);
          workingTargetPtr +=
            (stringSize > spaceLeft) ? spaceLeft : stringSize;
          spaceLeft -= (stringSize > spaceLeft) ? spaceLeft : stringSize;
          ++workingFormatPtr;
          break;
        default:
          va_end(varArgsList);
          std::sprintf(
            workingBuffer,
            "Exception(l_snprintf()): Format string %%%c not recognized.",
            *workingFormatPtr);
          throw std::invalid_argument(workingBuffer);
        }
      } else {
        *workingTargetPtr++ = *workingFormatPtr++;
        --spaceLeft;
      }
    }
    if(targetSize > 0) {
      *workingTargetPtr = '\0';
    }
    va_end(varArgsList);
    return static_cast<int>(workingTargetPtr - targetPtr);
  }

} // Anonymous namespace.


namespace brick {

  namespace common {

    // This constructor sets the internal "what()" message.  The
    Exception::
    Exception(const char* message)
      throw()
      : std::exception(),
        m_message(),
        m_payload(),
        m_payloadSize(0)
    {
      l_snprintf(m_message, BRICK_EXCEPTION_MESSAGE_LENGTH, "Exception: %s",
                 message);
    }


    // Constructor which accepts a few more details.
    Exception::
    Exception(const char* message, const char* fileName, int lineNumber)
      throw()
      : std::exception(),
        m_message(),
        m_payload(),
        m_payloadSize(0)
    {
      l_snprintf(m_message, BRICK_EXCEPTION_MESSAGE_LENGTH,
                 "Exception(%s, %d): %s",
                 fileName, lineNumber, message);
    }


    // Constructor which accepts a even more details.
    Exception::
    Exception(const char* message, const char* functionName,
              const char* fileName, int lineNumber)
      throw()
      : std::exception(),
        m_message(),
        m_payload(),
        m_payloadSize(0)
    {
      l_snprintf(m_message, BRICK_EXCEPTION_MESSAGE_LENGTH,
                 "Exception(%s, %s, %d): %s",
                 functionName, fileName, lineNumber, message);
    }


    // Copy constructor.
    Exception::
    Exception(const Exception& source)
      throw()
      : std::exception(),
        m_message(),
        m_payload(),
        m_payloadSize(0)
    {
      strncpy(m_message, source.m_message, BRICK_EXCEPTION_MESSAGE_LENGTH);
      m_message[BRICK_EXCEPTION_MESSAGE_LENGTH - 1] = '\0';
    }


    // This public method copies user-supplied data out of the
    // exception class.
    void
    Exception::
    getPayload(char* buffer, unsigned int& payloadSize) const throw()
    {
      payloadSize = m_payloadSize;
      memcpy(buffer, m_payload, m_payloadSize);
    }


    // This public method copies user-supplied data into the
    // exception class.
    void
    Exception::
    setPayload(char const* buffer, unsigned int bufferSize) throw()
    {
      this->setPayload(0, buffer, bufferSize);
    }


    // This public method copies user-supplied data into the
    // exception class.
    void
    Exception::
    setPayload(unsigned int skipBytes, char const* buffer,
               unsigned int bufferSize) throw()
    {
      skipBytes = std::min(
        skipBytes, static_cast<unsigned int>(BRICK_EXCEPTION_PAYLOAD_SIZE));
      if(skipBytes > m_payloadSize) {
        std::memset(m_payload + m_payloadSize, '\0', skipBytes - m_payloadSize);
      }
      bufferSize = std::min(
        bufferSize, BRICK_EXCEPTION_PAYLOAD_SIZE - skipBytes);
      memcpy(m_payload + skipBytes, buffer, bufferSize);
      m_payloadSize = skipBytes + bufferSize;
    }
      

    // Assignment operator.
    Exception& Exception::
    operator=(const Exception& source)
      throw()
    {
      if(this != &source) {
        strncpy(m_message, source.m_message, BRICK_EXCEPTION_MESSAGE_LENGTH);
        m_message[BRICK_EXCEPTION_MESSAGE_LENGTH - 1] = '\0';
        m_payloadSize = source.m_payloadSize;
        memcpy(m_payload, source.m_payload, m_payloadSize);
      }
      return *this;
    }


    // Protected constructor.
    Exception::
    Exception(const char* message, const char* childClassName)
      throw()
      : std::exception(),
        m_message(),
        m_payload(),
        m_payloadSize(0)
    {
      l_snprintf(m_message, BRICK_EXCEPTION_MESSAGE_LENGTH, "%s: %s",
                 childClassName, message);
    }


    // Protected constructor.
    Exception::
    Exception(const char* message, const char* childClassName,
              const char* functionName, const char* fileName,
              int lineNumber)
      throw()
      : std::exception(),
        m_message(),
        m_payload(),
        m_payloadSize(0)
    {
      if(functionName == 0) {
        l_snprintf(m_message, BRICK_EXCEPTION_MESSAGE_LENGTH,
                   "%s(%s, %d): %s",
                   childClassName, fileName, lineNumber, message);
      } else {
        l_snprintf(m_message, BRICK_EXCEPTION_MESSAGE_LENGTH,
                   "%s(%s, %s, %d): %s",
                   childClassName, functionName, fileName, lineNumber,
                   message);
      }
    }

  } // namespace common
  
} // namespace brick
