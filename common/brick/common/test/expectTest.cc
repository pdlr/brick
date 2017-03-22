/**
***************************************************************************
* @file inputStreamTest.cc
* 
* Source file defining tests for exception trace code.
*
* Copyright (C) 2005 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <sstream>
#include <string>
#include <brick/common/expect.hh>

namespace brick {

  namespace common {

    // We don't want to introduce a dependency on non-brick code for
    // unit testing, and the brick::test library is not available in
    // this context, so we just hack up some test functions.

    bool
    testExpect()
    {
      std::string inString;
      std::istringstream message0("Foo bar \t\n  baz  hork\n\nspork");
      message0 >> Expect("Foo") >> Expect(" bar ") >> Expect("\t\n  b")
               >> inString;
      if(!message0) {
        return false;
      }
      if(inString != "az") {
        return false;
      }
      message0 >> Expect("  ") >> Expect("hork\n\nspork");
      if(!message0) {
        return false;
      }

      std::istringstream message1("Foo bar \t\n  baz  hork\n\nspork");
      message1 >> Expect("Foo") >> Expect("bar");
      if(message1) {
        return false;
      }

      std::istringstream message2("Foo bar \t\n  baz  hork\n\nspork");
      message2 >> Expect("Foob");
      if(message2) {
        return false;
      }

      return true;
    }


    bool
    testExpectSkipWhitespace()
    {
      std::string inString;
      std::istringstream message0("Foo bar \t\n  baz  hork\n\nspork");
      Expect::FormatFlag flags(Expect::SkipWhitespace());
      message0 >> Expect("Foo", flags) >> Expect("bar ", flags)
               >> Expect("b", flags) >> inString;
      if(!message0) {
        return false;
      }
      if(inString != "az") {
        return false;
      }
      message0 >> Expect("hork", flags) >> Expect("spork", flags);
      if(!message0) {
        return false;
      }

      std::istringstream message1("Foo bar \t\n  baz  hork\n\nspork");
      message1 >> Expect("Foo", flags) >> Expect(" bar", flags);
      if(message1) {
        return false;
      }

      std::istringstream message2("Foo bar \t\n  baz  hork\n\nspork");
      message2 >> Expect("Foob", flags);
      if(message2) {
        return false;
      }

      return true;
    }


    bool
    testExpectSloppy()
    {
      std::string inString;
      std::istringstream message0("Foo bar \t\n  baz  hork\n\nspork");
      Expect::FormatFlag flags(Expect::Sloppy());
      message0 >> Expect("Foo", Expect::Sloppy()) >> Expect("nbyrn", flags)
               >> Expect("\t\n  b", flags) >> inString;
      if(!message0) {
        return false;
      }
      if(inString != "az") {
        return false;
      }
      message0 >> Expect("  ", flags) >> Expect("hork\n\nspork", flags);
      if(!message0) {
        return false;
      }

      std::istringstream message1("Foo bar \t\n  baz  hork\n\nspork");
      message1 >> Expect("Foo", flags) >> Expect(" bar", flags)
               >> Expect("123456", flags) >> inString;
      if(!message1) {
        return false;
      }
      if(inString != "az") {
        return false;
      }
      return true;
    }


    bool
    testExpectSloppySkipWhitespace()
    {
      std::string inString;
      std::istringstream message0("Foo bar \t\n  baz  hork\n\nspork");
      Expect::FormatFlag flags(Expect::SkipWhitespace() | Expect::Sloppy());
      message0 >> Expect("Foo", flags) >> Expect("nedl", flags)
               >> Expect("g", flags) >> inString;
      if(!message0) {
        return false;
      }
      if(inString != "az") {
        return false;
      }
      return true;
    }
    
  } // namespace common

} // namespace brick


// int main(int argc, char** argv)
int main(int, char**)
{
  bool result = true;
  result &= brick::common::testExpect();
  result &= brick::common::testExpectSkipWhitespace();
  result &= brick::common::testExpectSloppy();
  result &= brick::common::testExpectSloppySkipWhitespace();
  return (result ? 0 : 1);
}
