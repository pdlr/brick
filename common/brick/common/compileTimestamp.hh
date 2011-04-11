/**
***************************************************************************
* @file brick/common/compileTimestamp.hh
*
* Header file declaring CompileTimestamp class.
*
* Copyright (C) 2006-2010 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_COMPILETIMESTAMP_HH
#define BRICK_COMMON_COMPILETIMESTAMP_HH

namespace brick {

  namespace common {

    /**
     ** The CompileTimestamp class permits user code to conveniently
     ** assess when it was compiled.  To use this class, simply
     ** include compileTimestamp.hh and instantiate an instance of
     ** CompileTimestamp.  The resulting CompileTimestamp instance
     ** will reflect the time at which the source file was compiled.
     **/
    class CompileTimestamp
    {
    public:
      /** 
       * The default constructor sets the internal state to reflect the
       * time at which the compileTimestamp.hh> header was compiled.
       */
      CompileTimestamp()
        : m_day(0), m_hour(0), m_minute(0), m_month(0), m_second(0),
          m_year(0) {
        this->parseCompilerDateString(std::string(__DATE__));
        this->parseCompilerTimeString(std::string(__TIME__));
      }


      /** 
       * Destroys the CompileTimestamp instance.
       */
      ~CompileTimestamp() {}


      /** 
       * This member function returns a string of the form
       * "YYYY-MM-DD hh:mm:ss" indicating the local time at which the
       * file was compiled.
       * 
       * @return The return value is a string in ISO date/time format.
       */
      std::string
      getISOString();
    

    private:

      void
      parseCompilerDateString(const std::string& compilerTimeString);

    
      void
      parseCompilerTimeString(const std::string& compilerTimeString);

    
      int m_day;
      int m_hour;
      int m_minute;
      int m_month;
      int m_second;
      int m_year;
    };

  } // namespace common
    
} // namespace brick


#endif /* #ifndef BRICK_COMMON_COMPILETIMESTAMP_HH */
