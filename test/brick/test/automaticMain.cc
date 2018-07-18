/**
***************************************************************************
* @file brick/test/automaticMain.cpp
*
* Source file defining a main() routine to be used in unit testing.
*
* Copyright (C) 2007-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <set>
#include <brick/test/autoregister.hh>
#include <brick/test/runnableObject.hh>

namespace brick {

  namespace test {

    /// @cond privateCode
    namespace privateCode {

      // defined in autoregister.cpp
      std::set<RunnableObject*>&
      getGlobalTestSet();

    }
    /// @endcond

  }
}

using namespace brick::test;

int main(int, char**)
{
  bool aggregateResult = true;
  typedef std::set<RunnableObject*>::iterator TestFixtureIterator;
  std::set<RunnableObject*>& globalTestSet = privateCode::getGlobalTestSet();
  for(TestFixtureIterator iter = globalTestSet.begin();
      iter != globalTestSet.end();
      ++iter) {
    bool testResult = (*iter)->run();
    aggregateResult &= testResult;
  }
  return (aggregateResult ? 0 : 1);
}
