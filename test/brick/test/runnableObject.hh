/**
***************************************************************************
* @file brick/test/runnableObject.hh
*
* Header file declaring RunnableObject class, which is a parent class
* for the TestFixture class template.
*
* Copyright (C) 2007-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_TEST_RUNNABLEOBJECT_HH
#define BRICK_TEST_RUNNABLEOBJECT_HH

#include <brick/test/autoregister.hh>

namespace brick {

  namespace test {

    /**
     ** This class serves a parent class for the various
     ** TestFixture<foo> types, allowing virtual function dispatch to
     ** TestFixture<foo>::run(), and (if BRICK_TEST_USE_AUTOMATIC_MAIN
     ** is defined) automatic registration with the pre-written main()
     ** function.  User code will most likely never need to access
     ** this class.
     **/
    class RunnableObject {
    public:
      /**
       * Constructor registers with the pre-written main() function so
       * that it will be run automatically if the user chooses to link
       * with libbrickTestAutoMain.
       */
      RunnableObject() {registerTestFixture(*this);}


      /**
       * Destructor un-registers with the pre-written main() function
       * so that the no-longer-valid RunnableObject instance won't be
       * invoked automatically if the user chooses to link with
       * libbrickTestAutoMain.
       */
      virtual
      ~RunnableObject() {unregisterTestFixture(*this);}


      /**
       * Pure virtual run() method will dispatch to the run() method
       * of the subclass.
       *
       * @return The return value should indicate success or failure
       * of the test.
       */
      virtual bool
      run() = 0;

    };

  } // namespace test

} // namespace brick

#endif /* #ifdef BRICK_TEST_RUNNABLEOBJECT_HH */
