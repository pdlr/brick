/**
***************************************************************************
* @file brick/utilities/test/dateTest.cc
*
* Source file defining tests for date handling.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <sstream>
#include <string>

#include <brick/utilities/date.hh>
#include <brick/test/testFixture.hh>

using brick::test::TestFixture;

namespace brick {

  namespace utilities {

    class DateTest : public TestFixture<DateTest> {

    public:

      DateTest();
      ~DateTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__string();
      void testConstructor__size_t__size_t__size_t();
      void testConstructor__Date();
      void testDestructor();
      void testGetDay();
      void testGetMonth();
      void testGetYear();
      void testSetDay();
      void testSetMonth();
      void testSetYear();

      // Tests of non-member functions.
      void testOperatorEqualTo();
      void testOperatorGreaterThan();
      void testOperatorGreaterOrEqualTo();
      void testOperatorLessThan();
      void testOperatorLessOrEqualTo();
      void testOperatorNotEqualTo();
      void testOutputOperator();
      void testInputOperator();

    }; // class DateTest


    /* ============== Member Function Definititions ============== */

    DateTest::
    DateTest()
      : TestFixture<DateTest>("DateTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__string);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Date);
      BRICK_TEST_REGISTER_MEMBER(testDestructor);
      BRICK_TEST_REGISTER_MEMBER(testGetDay);
      BRICK_TEST_REGISTER_MEMBER(testGetMonth);
      BRICK_TEST_REGISTER_MEMBER(testGetYear);
      BRICK_TEST_REGISTER_MEMBER(testSetDay);
      BRICK_TEST_REGISTER_MEMBER(testSetMonth);
      BRICK_TEST_REGISTER_MEMBER(testSetYear);

      // Tests of non-member functions.
      BRICK_TEST_REGISTER_MEMBER(testOperatorEqualTo);
      BRICK_TEST_REGISTER_MEMBER(testOperatorGreaterThan);
      BRICK_TEST_REGISTER_MEMBER(testOperatorGreaterOrEqualTo);
      BRICK_TEST_REGISTER_MEMBER(testOperatorLessThan);
      BRICK_TEST_REGISTER_MEMBER(testOperatorLessOrEqualTo);
      BRICK_TEST_REGISTER_MEMBER(testOperatorNotEqualTo);
      BRICK_TEST_REGISTER_MEMBER(testOutputOperator);
      BRICK_TEST_REGISTER_MEMBER(testInputOperator);
    }


    void
    DateTest::
    testConstructor__void()
    {
      Date date0;
      BRICK_TEST_ASSERT(date0.getDay() == 0);
      BRICK_TEST_ASSERT(date0.getMonth() == 0);
      BRICK_TEST_ASSERT(date0.getYear() == 0);
    }

  
    void
    DateTest::
    testConstructor__string()
    {
      // Try a well-behaved input.
      Date date0("2003-1-11");
      BRICK_TEST_ASSERT(date0.getDay() == 11);
      BRICK_TEST_ASSERT(date0.getMonth() == 1);
      BRICK_TEST_ASSERT(date0.getYear() == 2003);

      // Try a not-so-well-behaved input, which (counting the extra
      // days), happens to be the same day.
      Date date1("2002-11-72");
      BRICK_TEST_ASSERT(date1.getDay() == date0.getDay());
      BRICK_TEST_ASSERT(date1.getMonth() == date0.getMonth());
      BRICK_TEST_ASSERT(date1.getYear() == date0.getYear());
    }

  
    void
    DateTest::
    testConstructor__size_t__size_t__size_t()
    {
      // Try a well-behaved input.
      Date date0(2002, 3, 01);
      BRICK_TEST_ASSERT(date0.getDay() == 1);
      BRICK_TEST_ASSERT(date0.getMonth() == 3);
      BRICK_TEST_ASSERT(date0.getYear() == 2002);

      // Try a not-so-well-behaved input which refers to the same day.
      // Note that 2002 is _not_ a leap year.
      Date date1(2002, 2, 29);
      BRICK_TEST_ASSERT(date1.getDay() == date0.getDay());
      BRICK_TEST_ASSERT(date1.getMonth() == date0.getMonth());
      BRICK_TEST_ASSERT(date1.getYear() == date0.getYear());

      // Make sure leap-years get noticed (2000) is a leap year.
      Date date2(2000, 2, 29);
      BRICK_TEST_ASSERT(date2.getDay() == 29);
      BRICK_TEST_ASSERT(date2.getMonth() == 2);
      BRICK_TEST_ASSERT(date2.getYear() == 2000);
    }

  
    void
    DateTest::
    testConstructor__Date()
    {
      Date date0("2003-1-11");
      Date date1(date0);
      BRICK_TEST_ASSERT(date1.getDay() == date0.getDay());
      BRICK_TEST_ASSERT(date1.getMonth() == date0.getMonth());
      BRICK_TEST_ASSERT(date1.getYear() == date0.getYear());
    }

  
    void
    DateTest::
    testDestructor()
    {
      // No explicit test.
    }

  
    void
    DateTest::
    testGetDay()
    {
      // No explicit test.
    }

  
    void
    DateTest::
    testGetMonth()
    {
      // No explicit test.
    }

  
    void
    DateTest::
    testGetYear()
    {
      // No explicit test.
    }

  
    void
    DateTest::
    testSetDay()
    {
      // Straightforward test.  No need to test any wacky day values,
      // since we've already verified that sanitize() is working fine.
      Date date0(2001, 1, 1);
      date0.setDay(12);
      BRICK_TEST_ASSERT(date0.getDay() == 12);
      BRICK_TEST_ASSERT(date0.getMonth() == 1);
      BRICK_TEST_ASSERT(date0.getYear() == 2001);
    }

  
    void
    DateTest::
    testSetMonth()
    {
      // Straightforward test.  No need to test any wacky day values,
      // since we've already verified that sanitize() is working fine.
      Date date0(2001, 1, 1);
      date0.setMonth(12);
      BRICK_TEST_ASSERT(date0.getDay() == 1);
      BRICK_TEST_ASSERT(date0.getMonth() == 12);
      BRICK_TEST_ASSERT(date0.getYear() == 2001);
    }

  
    void
    DateTest::
    testSetYear()
    {
      // Straightforward test.
      Date date0(2001, 1, 1);
      date0.setYear(2002);
      BRICK_TEST_ASSERT(date0.getDay() == 1);
      BRICK_TEST_ASSERT(date0.getMonth() == 1);
      BRICK_TEST_ASSERT(date0.getYear() == 2002);
    }

  
    // Tests of non-member functions.
    void
    DateTest::
    testOperatorEqualTo()
    {
      Date date0(2001, 2, 2);
      Date date1(2001, 2, 2);
      Date date2(2001, 2, 3);
      Date date3(2001, 3, 2);
      Date date4(2002, 2, 2);

      BRICK_TEST_ASSERT(date0 == date0);
      BRICK_TEST_ASSERT(date0 == date1);
      BRICK_TEST_ASSERT(!(date0 == date2));
      BRICK_TEST_ASSERT(!(date0 == date3));
      BRICK_TEST_ASSERT(!(date0 == date4));
      BRICK_TEST_ASSERT(!(date2 == date0));
      BRICK_TEST_ASSERT(!(date3 == date0));
      BRICK_TEST_ASSERT(!(date4 == date0));
    }

  
    void
    DateTest::
    testOperatorGreaterThan()
    {
      Date date0(2001, 2, 2);
      Date date1(2001, 2, 2);
      Date date2(2001, 2, 3);
      Date date3(2001, 3, 2);
      Date date4(2002, 2, 2);

      BRICK_TEST_ASSERT(!(date0 > date0));
      BRICK_TEST_ASSERT(!(date0 > date1));
      BRICK_TEST_ASSERT(!(date0 > date2));
      BRICK_TEST_ASSERT(!(date0 > date3));
      BRICK_TEST_ASSERT(!(date0 > date4));
      BRICK_TEST_ASSERT(date2 > date0);
      BRICK_TEST_ASSERT(date3 > date0);
      BRICK_TEST_ASSERT(date4 > date0);
    }

  
    void
    DateTest::
    testOperatorGreaterOrEqualTo()
    {
      Date date0(2001, 2, 2);
      Date date1(2001, 2, 2);
      Date date2(2001, 2, 3);
      Date date3(2001, 3, 2);
      Date date4(2002, 2, 2);

      BRICK_TEST_ASSERT(date0 >= date0);
      BRICK_TEST_ASSERT(date0 >= date1);
      BRICK_TEST_ASSERT(!(date0 >= date2));
      BRICK_TEST_ASSERT(!(date0 >= date3));
      BRICK_TEST_ASSERT(!(date0 >= date4));
      BRICK_TEST_ASSERT(date2 >= date0);
      BRICK_TEST_ASSERT(date3 >= date0);
      BRICK_TEST_ASSERT(date4 >= date0);
    }

  
    void
    DateTest::
    testOperatorLessThan()
    {
      Date date0(2001, 2, 2);
      Date date1(2001, 2, 2);
      Date date2(2001, 2, 3);
      Date date3(2001, 3, 2);
      Date date4(2002, 2, 2);

      BRICK_TEST_ASSERT(!(date0 < date0));
      BRICK_TEST_ASSERT(!(date0 < date1));
      BRICK_TEST_ASSERT(date0 < date2);
      BRICK_TEST_ASSERT(date0 < date3);
      BRICK_TEST_ASSERT(date0 < date4);
      BRICK_TEST_ASSERT(!(date2 < date0));
      BRICK_TEST_ASSERT(!(date3 < date0));
      BRICK_TEST_ASSERT(!(date4 < date0));
    }

  
    void
    DateTest::
    testOperatorLessOrEqualTo()
    {
      Date date0(2001, 2, 2);
      Date date1(2001, 2, 2);
      Date date2(2001, 2, 3);
      Date date3(2001, 3, 2);
      Date date4(2002, 2, 2);

      BRICK_TEST_ASSERT(date0 <= date0);
      BRICK_TEST_ASSERT(date0 <= date1);
      BRICK_TEST_ASSERT(date0 <= date2);
      BRICK_TEST_ASSERT(date0 <= date3);
      BRICK_TEST_ASSERT(date0 <= date4);
      BRICK_TEST_ASSERT(!(date2 <= date0));
      BRICK_TEST_ASSERT(!(date3 <= date0));
      BRICK_TEST_ASSERT(!(date4 <= date0));
    }

  
    void
    DateTest::
    testOperatorNotEqualTo()
    {
      Date date0(2001, 2, 2);
      Date date1(2001, 2, 2);
      Date date2(2001, 2, 3);
      Date date3(2001, 3, 2);
      Date date4(2002, 2, 2);

      BRICK_TEST_ASSERT(!(date0 != date0));
      BRICK_TEST_ASSERT(!(date0 != date1));
      BRICK_TEST_ASSERT(date0 != date2);
      BRICK_TEST_ASSERT(date0 != date3);
      BRICK_TEST_ASSERT(date0 != date4);
      BRICK_TEST_ASSERT(date2 != date0);
      BRICK_TEST_ASSERT(date3 != date0);
      BRICK_TEST_ASSERT(date4 != date0);
    }

  
    void
    DateTest::
    testOutputOperator()
    {
      Date date0(2001, 2, 2);
      Date date1;

      std::ostringstream outStream;
      outStream << date0;
      std::istringstream inStream(outStream.str());
      inStream >> date1;

      BRICK_TEST_ASSERT(date0 == date1);
    }

  
    void
    DateTest::
    testInputOperator()
    {
      // No explicit test.
    }

  } // namespace utilities
  
} // namespace brick


#ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION

int main(int argc, char** argv)
{
  brick::utilities::DateTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */

namespace {

  brick::utilities::DateTest currentTest;

}

#endif /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */
