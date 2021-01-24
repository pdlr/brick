/**
***************************************************************************
* @file brick/portability/timeUtilities.hh
*
* Header file declaring some useful time-related routines.
*
* (C) Copyright 2004-2011 David LaRose, dlr@davidlarose.com
* See accompanying LICENSE.TXT file for details.
*
***************************************************************************
**/

#ifndef BRICK_PORTABILITY_TIMEUTILITIES_HH
#define BRICK_PORTABILITY_TIMEUTILITIES_HH

#include <string>

namespace brick {

  namespace portability {

    enum TimeZone {
      BRICK_TZ_GMT,
      BRICK_TZ_LOCAL
    };


    /**
     * This function returns the current time as a double, the number of
     * seconds since some significant (and OS dependent) event in the
     * distant past.
     *
     * @return The return value is a double indicating the current time.
     */
    double
    getCurrentTime();


    /**
     * This function returns the current time as a string in the
     * format "YYYY-MM-DD HH:MM:SS.ss".
     *
     * @param tz This argument
     *
     * @return The return value
     */
    std::string
    getISOTimeString(TimeZone tz = BRICK_TZ_LOCAL);


    /**
     * This function causes the program to suspend execution for at
     * least as many seconds as specified by its argument.  The
     * resolution of the timer is implementation dependent, but will
     * always be some small fraction of a second.
     *
     * @param seconds This argument indicates how many seconds to
     * sleep.
     */
    void
    portableSleep(double seconds);

  } // namespace portability

} // namespace brick

#endif // #ifndef BRICK_PORTABILITY_TIMEUTILITIES_HH
