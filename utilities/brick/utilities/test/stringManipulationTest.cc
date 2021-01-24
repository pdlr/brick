/**
***************************************************************************
* @file brick/utilities/test/stringManipulationTest.cc
*
* Source file defining tests for string manipulation routines.
*
* Copyright (C) 2004-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <string>

#include <brick/utilities/stringManipulation.hh>
#include <brick/test/testFixture.hh>

using brick::test::TestFixture;

namespace brick {

  namespace utilities {

    class StringManipulationTest : public TestFixture<StringManipulationTest> {

    public:

      StringManipulationTest();
      ~StringManipulationTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Tests of member functions.
      void testCleanString();
      void testConvertString();
      void testJoinString();
      void testReplaceString();
      void testSplitString();
      void testStripString();
      void testWrapString();

    }; // class StringManipulationTest


    /* ============== Member Function Definititions ============== */

    StringManipulationTest::
    StringManipulationTest()
      : TestFixture<StringManipulationTest>("StringManipulationTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testCleanString);
      BRICK_TEST_REGISTER_MEMBER(testConvertString);
      BRICK_TEST_REGISTER_MEMBER(testJoinString);
      BRICK_TEST_REGISTER_MEMBER(testReplaceString);
      BRICK_TEST_REGISTER_MEMBER(testSplitString);
      BRICK_TEST_REGISTER_MEMBER(testStripString);
      BRICK_TEST_REGISTER_MEMBER(testWrapString);
    }


    void
    StringManipulationTest::
    testCleanString()
    {
      // The following string evaluates to "~foo\;bar\b;az".
      std::string testString = "~foo\\;bar\\b;az";

      // The following string evaluates to "\~foo\\\;bar\\b\;az".
      std::string resultString0 = "\\~foo\\\\\\;bar\\\\b\\;az";

      // The following string evaluates to "\~foo\;bar\b\;az".
      std::string resultString1 = "\\~foo\\;bar\\b\\;az";

      // Use cleanString to escape all instances of '\', ';', and '~'
      // using the quote character '\'.
      BRICK_TEST_ASSERT(cleanString(testString, "\\;~", '\\', false) ==
                        resultString0);

      // Use cleanString to escape all instances of '\', ';', and '~'
      // using the quote character '\', but don't quote things that are
      // already quoted.
      BRICK_TEST_ASSERT(cleanString(testString, "\\;~", '\\', true) ==
                        resultString1);
    }


    void
    StringManipulationTest::
    testConvertString()
    {
      // It would be nice if there were an easy way to query the type info
      // of built-in types.
      BRICK_TEST_ASSERT(sizeof(convertString<int>("11")) == sizeof(int));
      BRICK_TEST_ASSERT(sizeof(convertString<double>("11")) == sizeof(double));

      // We can, of course, check the values directly.
      BRICK_TEST_ASSERT(convertString<int>("11") == 11);
      BRICK_TEST_ASSERT(convertString<double>("11.5") == 11.5);
    }


    void
    StringManipulationTest::
    testJoinString()
    {
      std::vector<std::string> inputStringVector;
      std::string blankSeparator("");
      std::string unblankSeparator("HH");
      std::string string0("String0");
      std::string string1("String1");
      std::string string2("String2");
      std::string string3("String3");

      std::string outputString;
      outputString = joinString(inputStringVector);
      BRICK_TEST_ASSERT(outputString == "");
      outputString = joinString(inputStringVector, blankSeparator);
      BRICK_TEST_ASSERT(outputString == "");
      outputString = joinString(inputStringVector, unblankSeparator);
      BRICK_TEST_ASSERT(outputString == "");

      inputStringVector.push_back(string0);
      outputString = joinString(inputStringVector);
      BRICK_TEST_ASSERT(outputString == string0);
      outputString = joinString(inputStringVector, blankSeparator);
      BRICK_TEST_ASSERT(outputString == string0);
      outputString = joinString(inputStringVector, unblankSeparator);
      BRICK_TEST_ASSERT(outputString == string0);

      inputStringVector.push_back(string1);
      inputStringVector.push_back(string2);
      inputStringVector.push_back(string3);

      std::string referenceString;
      referenceString = string0 + string1 + string2 + string3;
      outputString = joinString(inputStringVector);
      BRICK_TEST_ASSERT(outputString == referenceString);
      outputString = joinString(inputStringVector, blankSeparator);
      BRICK_TEST_ASSERT(outputString == referenceString);

      referenceString = (string0 + unblankSeparator
                         + string1 + unblankSeparator
                         + string2 + unblankSeparator
                         + string3);
      outputString = joinString(inputStringVector, unblankSeparator);
      BRICK_TEST_ASSERT(outputString == referenceString);
    }


    void
    StringManipulationTest::
    testReplaceString()
    {
      std::string testString = "athisaisaaatestaastringaa";

      // Verify that simple replacement works OK.
      std::string newString = replaceString(testString, "a", "bbb");
      BRICK_TEST_ASSERT(newString
                        == "bbbthisbbbisbbbbbbbbbtestbbbbbbstringbbbbbb");

      // Verify that we're not faked out by overlapping instances of
      // target.
      newString = replaceString(testString, "aa", "bb");
      BRICK_TEST_ASSERT(newString == "athisaisbbatestbbstringbb");

      // One more check.
      newString = replaceString(testString, "aaa", "b");
      BRICK_TEST_ASSERT(newString == "athisaisbtestaastringaa");
    }


    void
    StringManipulationTest::
    testSplitString()
    {
      std::string testString = "athisaisaaatestaastringaa";

      // Test basic splitting.
      std::vector<std::string> stringParts = splitString(testString, "a");
      BRICK_TEST_ASSERT(stringParts.size() == 4);
      BRICK_TEST_ASSERT(stringParts[0] == "this");
      BRICK_TEST_ASSERT(stringParts[1] == "is");
      BRICK_TEST_ASSERT(stringParts[2] == "test");
      BRICK_TEST_ASSERT(stringParts[3] == "string");

      // Make sure "aaa" only counts as one "aa" occurance.
      stringParts = splitString(testString, "aa");
      BRICK_TEST_ASSERT(stringParts.size() == 3);
      BRICK_TEST_ASSERT(stringParts[0] == "athisais");
      BRICK_TEST_ASSERT(stringParts[1] == "atest");
      BRICK_TEST_ASSERT(stringParts[2] == "string");

      // Check again with a different separator.
      stringParts = splitString(testString, "aaa");
      BRICK_TEST_ASSERT(stringParts.size() == 2);
      BRICK_TEST_ASSERT(stringParts[0] == "athisais");
      BRICK_TEST_ASSERT(stringParts[1] == "testaastringaa");

      // Check that "includeNullStrings" works properly.
      stringParts = splitString(testString, "a", true);
      BRICK_TEST_ASSERT(stringParts.size() == 10);
      BRICK_TEST_ASSERT(stringParts[0] == "");
      BRICK_TEST_ASSERT(stringParts[1] == "this");
      BRICK_TEST_ASSERT(stringParts[2] == "is");
      BRICK_TEST_ASSERT(stringParts[3] == "");
      BRICK_TEST_ASSERT(stringParts[4] == "");
      BRICK_TEST_ASSERT(stringParts[5] == "test");
      BRICK_TEST_ASSERT(stringParts[6] == "");
      BRICK_TEST_ASSERT(stringParts[7] == "string");
      BRICK_TEST_ASSERT(stringParts[8] == "");
      BRICK_TEST_ASSERT(stringParts[9] == "");

      // Again to make sure that overlapping instances of separator
      // aren't a problem.
      stringParts = splitString(testString, "aa", true);
      BRICK_TEST_ASSERT(stringParts.size() == 4);
      BRICK_TEST_ASSERT(stringParts[0] == "athisais");
      BRICK_TEST_ASSERT(stringParts[1] == "atest");
      BRICK_TEST_ASSERT(stringParts[2] == "string");
      BRICK_TEST_ASSERT(stringParts[3] == "");

      // Make sure that "maxSplit" works as advertised.
      stringParts = splitString(testString, "a", false, 5);
      BRICK_TEST_ASSERT(stringParts.size() == 4);
      BRICK_TEST_ASSERT(stringParts[0] == "this");
      BRICK_TEST_ASSERT(stringParts[1] == "is");
      BRICK_TEST_ASSERT(stringParts[2] == "test");
      BRICK_TEST_ASSERT(stringParts[3] == "string");

      stringParts = splitString(testString, "a", false, 2);
      BRICK_TEST_ASSERT(stringParts.size() == 3);
      BRICK_TEST_ASSERT(stringParts[0] == "this");
      BRICK_TEST_ASSERT(stringParts[1] == "is");
      BRICK_TEST_ASSERT(stringParts[2] == "aatestaastringaa");

      // Make sure that "maxSplit" works as advertised when
      // "includeNullStrings" is specified.
      stringParts = splitString(testString, "a", true, 4);
      BRICK_TEST_ASSERT(stringParts.size() == 5);
      BRICK_TEST_ASSERT(stringParts[0] == "");
      BRICK_TEST_ASSERT(stringParts[1] == "this");
      BRICK_TEST_ASSERT(stringParts[2] == "is");
      BRICK_TEST_ASSERT(stringParts[3] == "");
      BRICK_TEST_ASSERT(stringParts[4] == "atestaastringaa");

      stringParts = splitString(testString, "a", true, 9);
      BRICK_TEST_ASSERT(stringParts.size() == 10);
      BRICK_TEST_ASSERT(stringParts[0] == "");
      BRICK_TEST_ASSERT(stringParts[1] == "this");
      BRICK_TEST_ASSERT(stringParts[2] == "is");
      BRICK_TEST_ASSERT(stringParts[3] == "");
      BRICK_TEST_ASSERT(stringParts[4] == "");
      BRICK_TEST_ASSERT(stringParts[5] == "test");
      BRICK_TEST_ASSERT(stringParts[6] == "");
      BRICK_TEST_ASSERT(stringParts[7] == "string");
      BRICK_TEST_ASSERT(stringParts[8] == "");
      BRICK_TEST_ASSERT(stringParts[9] == "");
    }


    void
    StringManipulationTest::
    testStripString()
    {
      std::string testString0 = "  \t \n this is\n\ttest string \t th\n \t\t  ";
      std::string testString1 = "this is\n\ttest string \t th";
      std::string testString2 = "is is\n\ttest string";

      // Test basic stripping.
      BRICK_TEST_ASSERT(stripString(testString0) == testString1);

      // Test stripping with user defined whitespace.
      BRICK_TEST_ASSERT(stripString(testString0, " \t\nth") == testString2);
    }


    void
    StringManipulationTest::
    testWrapString()
    {
      std::string testString0 =
        "This is a string to be wrapped.\nYesyesyes indeed it is";
      std::string referenceString0 =
        "This is\na string\nto be\nwrapped.\nYesyesyes\nindeed\nit is";
      std::string referenceString1 =
        "This is a\n> string\n> to be\n> wrapped.\n> Yesyesyes\n> indeed"
        "\n> it is";

      std::string testString2 =
        "    This is a short string to test an indexing bug.";

      std::string resultString0 =
        wrapString(testString0, "", 8, " \t\n", "\n");
      BRICK_TEST_ASSERT(resultString0 == referenceString0);

      std::string resultString1 =
        wrapString(testString0,  "> ", 10, " \t\n", "\n");
      BRICK_TEST_ASSERT(resultString1 == referenceString1);

      // Sencond test shouldn't change the string.
      std::string resultString2 = wrapString(testString2);
      BRICK_TEST_ASSERT(resultString2 == testString2);

    }

  } // namespace utilities

} // namespace brick


#ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION

int main(int argc, char** argv)
{
  brick::utilities::StringManipulationTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */

namespace {

  brick::utilities::StringManipulationTest currentTest;

}

#endif /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */
