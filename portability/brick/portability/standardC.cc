/**
***************************************************************************
* @file brick/portability/standardC.cpp
*
* Source file defining routines that conform to ANSI C, but are
* missing from various compilers.
*
* (C) Copyright 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying LICENSE file for details.
*
***************************************************************************
**/

#include <cstdarg>
#include <cstring>
#include <stdexcept>
#include <brick/portability/standardC.hh>

namespace brick {

  namespace portability {

#if HAVE_SNPRINTF
#else

    int snprintf(char* targetPtr, size_t targetSize, const char* formatPtr,
                 ...) {
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
              "Exception(brick::snprintf()): Format string %%%c not recognized.",
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

#endif /* #if HAVE_SNPRINTF */

  } // namespace portability

} // namespace brick
