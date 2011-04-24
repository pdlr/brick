/**
***************************************************************************
* @file brick/utilities/timeUtilities.hh
*
* Header file declaring some useful time-related routines.
*
* (C) Copyright 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying LICENSE file for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_TIMEUTILITIES_HH
#define BRICK_UTILITIES_TIMEUTILITIES_HH

#include <brick/portability/timeUtilities.hh>

namespace brick {

  namespace utilities {
    
    /** 
     * This function returns the current time as a double, the number of
     * seconds since some significant (and OS dependent) event in the
     * distant past.
     */
    using portability::getCurrentTime;

    
    /** 
     * This function causes the program to suspend execution for at
     * least as many seconds as specified by its argument.  The
     * resolution of the timer is implimentation dependent, will
     * always be some small fraction of a second.
     */
    using portability::portableSleep;


    /**
     ** This class represents a timer, useful for measuring time
     ** interals.
     **/
    class Timer {
    public:

      /** 
       * This constructor initializes the timer and starts timing.
       * 
       * @param timeUnit This argument 
       */
      Timer()
	: m_startTime(0.0) {this->reset();}

    
      /** 
       * Destructor.
       */
      ~Timer() {}


      /** 
       * The function calculates and returns the time since the last
       * reset (or since the timer was created).
       * 
       * @return The return value is the elapsed time in seconds.
       */
      double
      getElapsedTime() {return getCurrentTime() - m_startTime;}


      /** 
       * This function behaves just like getElapsedTime(), but also
       * resets the timer to zero.
       * 
       * @return The return value is the elapsed time in seconds.
       */
      double
      reset() {
        double t0 = m_startTime;
        m_startTime = getCurrentTime();
        return m_startTime - t0;
      }

    private:
      
      double m_startTime;
      
    }; // class Timer
    
  } // namespace utilities
  
} // namespace brick

#endif /* #ifndef BRICK_UTILITIES_TIMEUTILITIES_HH */
