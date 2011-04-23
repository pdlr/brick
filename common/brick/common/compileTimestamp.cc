/**
***************************************************************************
* @file brick/common/compileTimestamp.cc
*
* Source file defining CompileTimestamp class.
*
* Copyright (C) 2006-2010 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <iomanip>
#include <sstream>
#include <brick/common/compileTimestamp.hh>
#include <brick/common/exception.hh>
#include <brick/common/expect.hh>

namespace brick {

  namespace common {

    std::string
    CompileTimestamp::
    getISOString()
    {
      std::ostringstream isoStream;
      isoStream << std::setfill('0')
		<< std::setw(4) << m_year << "-"
		<< std::setw(2) << m_month + 1 << "-"
		<< std::setw(2) << m_day + 1 << " "
		<< std::setw(2) << m_hour << ":"
		<< std::setw(2) << m_minute << ":"
		<< std::setw(2) << m_second;
      return isoStream.str();
    }
    

    void
    CompileTimestamp::
    parseCompilerDateString(const std::string& compilerDateString)
    {
      std::istringstream compilerDateStream(compilerDateString);

      std::string monthString;
      compilerDateStream >> monthString;
      if(monthString == "Jan") {m_month = 0;}
      else if(monthString == "Feb") {m_month = 1;}
      else if(monthString == "Mar") {m_month = 2;}
      else if(monthString == "Apr") {m_month = 3;}
      else if(monthString == "May") {m_month = 4;}
      else if(monthString == "Jun") {m_month = 5;}
      else if(monthString == "Jul") {m_month = 6;}
      else if(monthString == "Aug") {m_month = 7;}
      else if(monthString == "Sep") {m_month = 8;}
      else if(monthString == "Oct") {m_month = 9;}
      else if(monthString == "Nov") {m_month = 10;}
      else if(monthString == "Dec") {m_month = 11;}
      else {
	std::ostringstream message;
	message << "Unable to recover month from date string \""
		<< compilerDateString << "\".";
	BRICK_THROW(ValueException,
		  "CompileTimestamp::parseCompilerDateString()",
		  message.str().c_str());
      }

      compilerDateStream >> m_day;
      --m_day;

      compilerDateStream >> m_year;

      if(!compilerDateStream) {
	std::ostringstream message;
	message << "Unable to parse date string \"" << compilerDateString
		<< "\".";
	BRICK_THROW(ValueException,
		  "CompileTimestamp::parseCompilerDateString()",
		  message.str().c_str());
      }
    }

  
    void
    CompileTimestamp::
    parseCompilerTimeString(const std::string& compilerTimeString)
    {
      std::istringstream inputStream(compilerTimeString);
      inputStream >> m_hour >> Expect(":")
                  >> m_minute >> Expect(":")
                  >> m_second;
      if(!inputStream) {
	std::ostringstream message;
	message << "Unable to parse time string \"" << compilerTimeString
		<< "\".";
	BRICK_THROW(ValueException,
		  "CompileTimestamp::parseCompilerTimeString()",
		  message.str().c_str());
      }
    }
  
  } // namespace common

} // namespace brick
