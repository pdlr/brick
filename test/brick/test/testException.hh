/**
***************************************************************************
* @file brick/test/testException.hh
*
* Header file declaring a specialization of brick::Exception which will
* be thrown to indicate failed unit tests.
*
* Copyright (C) 2003-2004, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_TEST_TESTEXCEPTION_HH
#define BRICK_TEST_TESTEXCEPTION_HH

#include <brick/common/exception.hh>


namespace brick {

  /**
   ** This namespace contains a unit testing library.
   **/
  namespace test {

    /**
     ** This exception is thrown to indicate a failed test.  It is
     ** usually thrown by one of the BRICK_TEST_ASSERT() macros, but
     ** it can also be thrown explicitly by the test code.
     **/
    class TestException;
    BRICK_DECLARE_EXCEPTION_TYPE(TestException, brick::common::Exception);

  } // namespace test

} // namespace brick

#endif /* #ifndef BRICK_TEST_TESTEXCEPTION_HH */
