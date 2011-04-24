/**
***************************************************************************
* @file brick/test/autoregister.cpp
*
* Source file allowing control of what gets run by the main() routine
* in libdlrTestAutomaticMain.
*
* Copyright (C) 2007-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_TEST_USE_AUTOMATIC_MAIN
#define BRICK_TEST_USE_AUTOMATIC_MAIN
#endif

#include <set>
#include <brick/test/autoregister.hh>
#include <brick/test/runnableObject.hh>

namespace brick {

  namespace test {

    /// @cond privateCode
    namespace privateCode {

      std::set<RunnableObject*>&
      getGlobalTestSet() {
	static std::set<RunnableObject*> globalTestSet;
	return globalTestSet;
      }
      
    }
    /// @endcond


    // Schedule a TestFixture to be run automatically.
    void
    private_registerTestFixture(RunnableObject& testFixture)
    {
      (privateCode::getGlobalTestSet()).insert(&testFixture);
    }


    // Un-schedule a TestFixture so that it will no longer be run
    // automatically.
    void
    private_unregisterTestFixture(RunnableObject& testFixture)
    {
      (privateCode::getGlobalTestSet()).erase(&testFixture);
    }

  }
}
