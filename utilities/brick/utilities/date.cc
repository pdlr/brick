/**
***************************************************************************
* @file brick/utilities/date.cc
*
* Source file defining Date class.
*
* Copyright (C) 2003-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <iostream>
#include <sstream>
#include <iomanip>
#include <brick/utilities/stringManipulation.hh>
#include <brick/utilities/date.hh>

namespace brick {

  namespace utilities {

    /// Set from an ISO date string
    void
    Date::
    setISODate(std::string const& isoString)
    {
      // Be a little forgiving of format.  Replace both "-" and "/"
      // with whitespace (not just "-")
      std::string string0 = replaceString(isoString, "/", " ");
      std::string string1 = replaceString(isoString, "-", " ");

      // Now parse string.
      std::istringstream inputStream(string1);
      inputStream >> m_year >> m_month >> m_day;
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't parse ISO date string: " << isoString << std::endl;
        BRICK_THROW(brick::common::ValueException,
                    "Date::setISODate()", message.str().c_str());
      }

      // Make sure month and day are sane.
      this->sanitize();
    }

    // Returns the number of days in a given month.
    size_t Date::
    getDayCount(size_t month, size_t year)
    {
      switch(month) {
      case 1:
      case 3:
      case 5:
      case 7:
      case 8:
      case 10:
      case 12:
        return 31;
        break;
      case 4:
      case 6:
      case 9:
      case 11:
        return 30;
        break;
      case 2:
        if(this->isLeapYear(year)) {
          return 29;
        }
        return 28;
        break;
      default:
        std::ostringstream message;
        message << "Invalid month: " << month << std::endl;
        BRICK_THROW(brick::common::ValueException, "Date::getDayCount()",
                    message.str().c_str());
        break;
      }
      BRICK_THROW(brick::common::LogicException, "Date::getDayCount()",
                  "Faulty switch logic\n");
      return 0;
    }

    // Chech for leapyear.
    bool Date::
    isLeapYear(size_t year)
    {
      if(year % 4 == 0) {
        return true;
      }
      return false;
    }

    // Make sure day and date numbers are calender-possible.
    void Date::
    sanitize()
    {
      while(1) {
        // Fix month first, since this might bump us into a leapyear
        while(m_month > 12) {
          m_month -= 12;
          m_year += 1;
        }
        // Done?
        if(m_day <= getDayCount(m_month, m_year)) {
          break;
        }
        m_day -= getDayCount(m_month, m_year);
        m_month++;
      }
    }

    bool
    operator<(Date const& date0, Date const& date1)
    {
      if(date0.getYear() > date1.getYear()) {
        return false;
      }
      if(date0.getYear() == date1.getYear()) {
        if(date0.getMonth() > date1.getMonth()) {
          return false;
        }
        if(date0.getMonth() == date1.getMonth()) {
          if(date0.getDay() >= date1.getDay()) {
            return false;
          }
        }
      }
      return true;
    }

    bool
    operator>(Date const& date0, Date const& date1)
    {
      return (!((date0 == date1) || (date0 < date1)));
    }

    bool
    operator<=(Date const& date0, Date const& date1)
    {
      return ((date0 == date1) || (date0 < date1));
    }

    bool
    operator>=(Date const& date0, Date const& date1)
    {
      return !(date0 < date1);
    }

    bool
    operator==(Date const& date0, Date const& date1)
    {
      return ((date0.getYear() == date1.getYear())
              && (date0.getMonth() == date1.getMonth())
              && (date0.getDay() == date1.getDay()));
    }

    bool
    operator!=(Date const& date0, Date const& date1)
    {
      return !(date0 == date1);
    }

    std::ostream&
    operator<<(std::ostream& stream, Date const& date)
    {
      stream << std::setfill('0') << std::setw(4) << date.getYear() << "-"
             << std::setw(2) << date.getMonth() << "-"
             << std::setw(2) << date.getDay() << std::setfill(' ');
      return stream;
    }

    std::istream&
    operator>>(std::istream& stream, Date& date)
    {
      std::string isoString;
      stream >> isoString;

      // Can't read from stream if it's in a bad state.
      if(!stream) {
        return stream;
      }

      try {
        Date newDate(isoString);	// Make sure things parse OK before
        date = newDate;		// assigning to argument date.
      } catch(const brick::common::ValueException&) {
        stream.clear(std::ios_base::failbit);
        return stream;
      }

      return stream;
    }


  } // namespace utilities

} // namespace brick
