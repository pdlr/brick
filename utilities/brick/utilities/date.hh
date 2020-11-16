/**
***************************************************************************
* @file brick/utilities/date.hh
*
* Header file declaring Date class.
*
* Copyright (C) 2003-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_DATE_HH
#define BRICK_UTILITIES_DATE_HH

namespace brick {

  namespace utilities {

    /**
     ** As you'd expect, the Date class represents a date, like 2003-03-15.
     ** Note that the format is yyyy-mm-dd.
     **/
    class Date {
    public:
      /**
       * The default constructor initializes to the imaginary date
       * 0000-00-00.
       */
      Date()
	: m_day(0), m_month(0), m_year(0) {}


      /**
       * This constructor sets the date from ISO string "yyyy-mm-dd".
       *
       * @param isoString This argument specifies the date.
       */
      Date(std::string const& isoString) {this->setISODate(isoString);}


      /**
       * This constructor sets the date explicitly.
       *
       * @param year This argument specifies the year.
       *
       * @param month This argument specifies the month.
       *
       * @param day This argument specifies the day.
       */
      Date(size_t year, size_t month, size_t day)
	: m_day(day), m_month(month), m_year(year) {this->sanitize();}


      /**
       * The copy constructor deep copies its argument.
       *
       * @param source This argument is the Date instance to be copied.
       */
      Date(Date const& source)
	: m_day(source.m_day), m_month(source.m_month),
	  m_year(source.m_year) {}


      /**
       * Destructor.
       */
      ~Date() {}


      /**
       * The copy assignment operator deep copies its argument.
       *
       * @param source This argument is the Date instance to be copied.
       */
      Date&
      operator=(Date const& source)
      {
          m_day = source.m_day;
          m_month = source.m_month;
          m_year = source.m_year;
          return *this;
      }


      /**
       * The function returns the day.
       *
       * @return The return value is the day.
       */
      size_t
      getDay() const {return m_day;}


      /**
       * The function returns the month.
       *
       * @return The return value is the month.
       */
      size_t
      getMonth() const {return m_month;}


      /**
       * The function returns the year.
       *
       * @return The return value is the year.
       */
      size_t
      getYear() const {return m_year;}


      /**
       * This function explicitly sets the day.
       *
       * @param day This argument is the specified day.
       */
      void
      setDay(size_t day) {m_day = day; this->sanitize();}


      /**
       * This function sets the date from ISO string "yyyy-mm-dd".
       *
       * @param isoString This argument specifies the date.
       */
      void
      setISODate(std::string const& isoString);


      /**
       * This function explicitly sets the month.
       *
       * @param month This argument is the specified month.
       */
      void
      setMonth(size_t month) {m_month = month; this->sanitize();}


      /**
       * This function explicitly sets the year.
       *
       * @param year This argument is the specified year.
       */
      void
      setYear(size_t year) {m_year = year;}

    private:
      // Private member functions.

      /// Returns the number of days in a given month.
      size_t
      getDayCount(size_t month, size_t year);

      // Check for leapyear.
      bool
      isLeapYear(size_t year);

      /// Make sure day and month numbers are calendar-possible.
      void
      sanitize();

      // Private data members.
      size_t m_day;
      size_t m_month;
      size_t m_year;
    }; // class Date


    /**
     * This function compares two dates.
     *
     * @param date0 This argument is the first date to compare.
     *
     * @param date1 This argument is the second date to compare.
     *
     * @return The return value is true if first date is earlier than
     * the second, otherwise false.
     */
    bool
    operator<(Date const& date0, Date const& date1);


    /**
     * This function compares two dates.
     *
     * @param date0 This argument is the first date to compare.
     *
     * @param date1 This argument is the second date to compare.
     *
     * @return The return value is true if first date is later than
     * the second, otherwise false.
     */
    bool
    operator>(Date const& date0, Date const& date1);


    /**
     * This function compares two dates.
     *
     * @param date0 This argument is the first date to compare.
     *
     * @param date1 This argument is the second date to compare.
     *
     * @return The return value is true if first date is earlier than or
     * equal to the second, otherwise false.
     */
    bool
    operator<=(Date const& date0, Date const& date1);


    /**
     * This function compares two dates.
     *
     * @param date0 This argument is the first date to compare.
     *
     * @param date1 This argument is the second date to compare.
     *
     * @return The return value is true if first date is later than or
     * equal to the second, otherwise false.
     */
    bool
    operator>=(Date const& date0, Date const& date1);


    /**
     * This function compares two dates.
     *
     * @param date0 This argument is the first date to compare.
     *
     * @param date1 This argument is the second date to compare.
     *
     * @return The return value is true if first date is the same as the
     * second, otherwise false.
     */
    bool
    operator==(Date const& date0, Date const& date1);


    /**
     * This function compares two dates.
     *
     * @param date0 This argument is the first date to compare.
     *
     * @param date1 This argument is the second date to compare.
     *
     * @return The return value is true if first date is different than
     * the second, otherwise false.
     */
    bool
    operator!=(Date const& date0, Date const& date1);


    /**
     * This operator formats a date in ISO standard format
     * ("yyyy-mm-dd") for stream output.
     *
     * @param stream This argument is the stream to which the format is
     * to be written.
     *
     * @param date This argument is the date to be formatted.
     *
     * @return The return value is the stream, after the output
     * operation.
     */
    std::ostream&
    operator<<(std::ostream& stream, Date const& date);


    /**
     * This operator reads a date in ISO standard format ("yyyy-mm-dd")
     * from an input stream.
     *
     * @param stream This argument is the stream from which the date is
     * to be read.
     *
     * @param date This argument is the date to be modified.
     *
     * @return The return value is the stream, after the input
     * operation.
     */
    std::istream&
    operator>>(std::istream& stream, Date& date);

  } // namespace utilities

} // namespace brick

#endif /* #ifndef BRICK_UTILITIES_DATE_HH */
