/**
***************************************************************************
* @file brick/test/autoregister.hh
*
* Header file allowing control of what gets run by the main() routine
* in libbrickTestAutomaticMain.
*
* Copyright (C) 2007-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_TEST_AUTOREGISTER_HH
#define BRICK_TEST_AUTOREGISTER_HH

namespace brick {

  namespace test {

    // Forward declaration.
    class RunnableObject;


    void
    private_registerTestFixture(RunnableObject& testFixture);

    void
    private_unregisterTestFixture(RunnableObject& testFixture);


#ifndef BRICK_TEST_NO_AUTOMATIC_REGISTRATION

    /**
     * This function Registers a test fixture with the pre-written
     * main() function so that it will be run automatically if the
     * user chooses to link with libbrickTestAutoMain.  This function
     * can be disabled by defining BRICK_TEST_NO_AUTOMATIC_REGISTRATION.
     *
     * @param testFixture This argument is the test fixture to be
     * registered.
     */
    inline void
    registerTestFixture(RunnableObject& testFixture) {
      private_registerTestFixture(testFixture);
    }


    /**
     * This function unregisters a test fixture with the pre-written
     * main() function so that it will not be be run automatically if
     * the user chooses to link with libbrickTestAutoMain.  It undoes
     * the effect of registerTestFixture().  This function can be
     * disabled by defining BRICK_TEST_NO_AUTOMATIC_REGISTRATION.
     *
     * @param testFixture This argument is the test fixture to be
     * registered.
     */
    inline void
    unregisterTestFixture(RunnableObject& testFixture) {
      private_unregisterTestFixture(testFixture);
    }

#else /* #ifndef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */

    inline void
    registerTestFixture(RunnableObject&) {}

    inline void
    unregisterTestFixture(RunnableObject&) {}

#endif /* #ifndef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */

  } // namespace test

} // namespace brick

#endif /* #ifdef BRICK_TEST_AUTOREGISTER_HH */
