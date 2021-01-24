/**
***************************************************************************
* @file brick/portability/timeUtilities.cpp
*
* Source file declaring some useful time-related routines.
*
* (C) Copyright 2004-2011 David LaRose, dlr@davidlarose.com
* See accompanying LICENSE.TXT file for details.
*
***************************************************************************
**/

#include <brick/common/exception.hh>
#include <brick/portability/timeUtilities.hh>

/* ===================== Common includes ===================== */


/* ===================== End common includes ===================== */

#ifdef _WIN32

/* ===================== Windows includes ===================== */

#include <time.h>
#include <windows.h>
#include <sys/timeb.h>

/* ===================== End Windows includes ===================== */

#else /* #ifdef _WIN32 */

/* ===================== Linux includes ===================== */

#include <time.h>
#include <sys/time.h>
#include <sys/select.h>

/* ===================== End Linux includes ===================== */

#endif /* #ifdef _WIN32 */


/* ===================== Common code ===================== */

namespace brick {

  namespace portability {

    // Pass.

  } // namespace portability

} // namespace brick

/* ===================== End common code ===================== */


#ifdef _WIN32

/* ===================== Windows code ===================== */

// Anonymous namespace for stuff local to this file.
namespace brick {

  namespace portability {

    double
    getCurrentTime() {
      struct __timeb64 tp;
      _ftime64(&tp);
      return tp.time + (tp.millitm * 0.001);
    }


    void
    portableSleep(double seconds) {
      if(seconds > 0.0) {
        int milliseconds = static_cast<int>(1000.0 * seconds + 0.5);
        Sleep(milliseconds);
      }
    }

  } // namespace portability

} // namespace brick

/* ===================== End Windows code ===================== */

#else /* #ifdef _WIN32 */

/* ===================== Linux code ===================== */

#include <cmath>

namespace brick {

  namespace portability {

    double
    getCurrentTime() {
      struct timeval tp;
      gettimeofday(&tp, 0);
      return tp.tv_sec + (tp.tv_usec * 0.000001);
    }


    std::string
    getISOTimeString(TimeZone timeZone) {
      // 20 characters is enough for "YYYY-MM-DD hh:mm:ss\0".
      const int bufferSize = 20;
      char buffer[bufferSize];
      time_t timeSinceEpoch;
      time_t returnValue = time(&timeSinceEpoch);
      if(returnValue == -1) {
        BRICK_THROW(brick::common::RunTimeException, "getISOTimeString()",
                    "Call to time() failed.");
      }
      struct tm brokenDownTime;
      switch(timeZone) {
      case BRICK_TZ_GMT:
        gmtime_r(&timeSinceEpoch, &brokenDownTime);
        break;
      case BRICK_TZ_LOCAL:
        tzset();
        localtime_r(&timeSinceEpoch, &brokenDownTime);
        break;
      default:
        BRICK_THROW(brick::common::LogicException, "getISOTimeString()",
                    "Unexpected value in switch() argument.");
        break;
      }

      strftime(buffer, bufferSize, "%Y-%m-%d %H:%M:%S", &brokenDownTime);
      return std::string(buffer);
    }


    void
    portableSleep(double seconds) {
      if(seconds > 0.0) {
        double integerPart;
        double fractionalPart = std::modf(seconds, &integerPart);
        struct timeval timeout;
        timeout.tv_sec = static_cast<int>(integerPart);
        timeout.tv_usec =
          static_cast<int>(1000000.0 * fractionalPart + 0.5);
        if(timeout.tv_usec == 1000000) {
          ++timeout.tv_sec;
          timeout.tv_usec = 0;
        }
        select(0, NULL, NULL, NULL, &timeout);
      }
    }

  } // namespace portability

} // namespace brick

/* ===================== End Linux code ===================== */

#endif /* #ifdef _WIN32 */
