/**
***************************************************************************
* @file brick/utilities/frequencyGoverner.hh
*
* Header file declaring the FrequencyGoverner class.
*
* Copyright (C) 2007-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_FREQUENCYGOVERNER_HH
#define BRICK_UTILITIES_FREQUENCYGOVERNER_HH

#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace utilities {

    /**
     ** The FrequencyGoverner class allows you to conveniently throttle
     ** a loop so that it runs at a specific speed.  Here's an example
     ** of how it might be used:
     **
     ** @code
     **   double frequency = 30.0; // Hz.
     **   double duration = 10.0;  // Seconds.
     **   FrequencyGoverner governer(frequency, duration);
     **   while(!governer.isFinished()) {
     **     // Do something.
     **     // [..]
     **     governer.sleepIfNecessary();
     **   }
     ** @endcode
     **/
    class FrequencyGoverner {
    public:

      /**
       * The constructor initializes the FrequencyGoverner and starts
       * the internal timer.
       *
       * @param frequency This argument specifies how fast you'd like
       * the loop to run.  Setting this argument to zero indicates no
       * throttling of loop speed.
       *
       * @param duration This argument specifies how long you'd like the
       * loop to run.  Setting this argument to zero indicates no
       * timeout.
       */
      explicit
      FrequencyGoverner(double frequency, double duration = 0.0)
        : m_count(0), m_duration(duration), m_period(1.0 / frequency),
          m_startTime(getCurrentTime()) {}


      /**
       * The destructor destroys the FrequencyGoverner instance and deletes any
       * file created by the constructor, thereby releasing the lock.
       */
      ~FrequencyGoverner() {}


      /**
       * This method returns the average frequency at which method
       * sleepIfNecessary() has been called.  This average includes
       * any time spent between the constructor and the first call to
       * sleepIfNecessary().  As long as the other code in the
       * frequency-governed loop isn't too heavy, the return value
       * should be approximately equal to the frequency specified in
       * the constructor call.
       *
       * @return The return value is an estimage of the actual
       * frequency at which the sleepIfNecessary is being called.
       */
      double
      getActualFrequency() {
        if(m_count == 0) {return 0.0;}
        return (getCurrentTime() - m_startTime) / m_count;
      }


      /**
       * This method returns the number of times that method
       * sleepIfNecessary() has been called.
       *
       * @return The return value indicates the totall number of calls
       * to sleepIfNecessary().
       */
      unsigned int
      getCount() {return static_cast<unsigned int>(m_count);}


      /**
       * This method returns how long, in seconds, method
       * sleepIfNecessary() would have slept if it had been called
       * instead of getSlack().
       *
       * @return The return value is positive if we are ahead of
       * schedule, negative otherwise.
       */
      double
      getSlack() {return m_startTime + m_count * m_period - getCurrentTime();}


      /**
       * This method reports whether or not the specified period (from
       * constructor argument "duration") has elapsed.
       *
       * @return The return value is true if the FrequencyGoverner was
       * instantiated more than "duration" seconds ago, false otherwise.
       */
      inline bool
      isFinished() {
        return ((m_duration != 0.0)
                && ((getCurrentTime() - m_startTime) > m_duration));
      }


      /**
       * This method sleeps just long enough to prevent a loop from
       * running faster than the specifed frequency.
       */
      inline void
      sleepIfNecessary() {
        // Function portableSleep is smart enough to return
        // immediately for negative arguments.
        portableSleep(this->getSlack());
        ++m_count;
      }

    private:

      int m_count;
      double m_duration;
      double m_period;
      double m_startTime;

    }; // class FrequencyGoverner

  } // namespace utilities

} // namespace brick

#endif /* #ifndef BRICK_UTILITES_FREQUENCYGOVERNER_HH */
