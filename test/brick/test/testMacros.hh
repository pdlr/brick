/**
***************************************************************************
* @file brick/test/testMacros.hh
* 
* Header file declaring macros for brick::test library.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_TEST_TESTMACROS_HH
#define BRICK_TEST_TESTMACROS_HH

// Anyone who is using these macros will also want to have included
// testException.h .
#include <brick/test/testException.hh>

/**
 ** This macro throws a TestException if its argument does not
 ** evaluate to true.
 **
 ** Example: BRICK_TEST_ASSERT(x == 10);
 **/
#define BRICK_TEST_ASSERT(assertion) { \
  if(!(assertion)) { \
    BRICK_THROW2(brick::test::TestException, (#assertion)); \
  } \
}


/**
 ** This macro throws a TestException if its argument does not throw
 ** the specified exception when evaluated.
 **
 ** Example: BRICK_TEST_ASSERT_EXCEPTION(ValueException, functionWhichThrows());
 **/
#define BRICK_TEST_ASSERT_EXCEPTION(exceptionType, assertion) { \
  try { \
    assertion; \
    BRICK_THROW2(brick::test::TestException, \
               (#assertion " should throw " #exceptionType)); \
  } catch(const exceptionType&) {} \
}


/** 
 ** This macro makes it easier to call the
 ** TestFixture<>::registerTest() member function.  It only makes
 ** sense to use this macro from within a member function of a
 ** subclass of TestFixture.
 **/
#define BRICK_TEST_REGISTER_MEMBER(test) { \
  this->registerTest((#test), &TestFixtureType::test); \
}


#endif /* #ifndef BRICK_TEST_TESTMACROS_HH */
