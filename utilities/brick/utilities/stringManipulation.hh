/**
***************************************************************************
* @file brick/utilities/stringManipulation.hh
*
* Header file declaring routines for working with strings.
*
* Copyright (C) 2003-2011, David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_STRINGMANIPULATION_HH
#define BRICK_UTILITIES_STRINGMANIPULATION_HH

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <brick/common/exception.hh>
#include <brick/utilities/exception.hh>

namespace brick {

  namespace utilities {

    /**
     * This function takes an input string and returns a string in which
     * all special shell characters have been escaped.  The list of
     * special characters is specified as a string, and defaults to the
     * following characters, including the double quote character:
     * "~#$&*()\|[]{};'`<>/?".  The third argument specifies the quote
     * character, and defaults to the backslash character ("\").  The
     * third argument affects how quote characters in the input are
     * handled.  By default, quote characters are escaped just like any
     * other special character.  If the third argument is set to true,
     * then it is assumed that parts of the string have already been
     * quoted, and quote characters are not quoted.  For example:
     *
     *   cleanString("~foo\;bar\b;az", "\;~", '\', false);
     *   will return the string "\~foo\\\;bar\\b\;az",
     *
     * and
     *
     *   cleanString("~foo\;bar\b;az", "\;~", '\', true);
     *   will return the string "\~foo\;bar\b\;az".
     *
     * @param inputString This is the string to be cleaned.
     *
     * @return The return value is a copy of inputString in which
     * special characters are quoted.
     */
    std::string
    cleanString(std::string const& inputString,
                std::string const&
                specialCharacters="\"~#$&*()\\|[]{};'`<>/?",
                char quoteCharacter='\\',
                bool alreadyQuoted=false);


    /**
     * This function converts the string inputString to Type.  For
     * example:
     *
     * @code
     *   int intValue = convertString<int>("23");
     * @endcode
     *
     * Conversion is done using the stream input operator of the
     * target class.  If the conversion fails, a
     * brick::utilities::ConversionException is thrown.
     *
     * @param inputString This argument is the string to be converted.
     *
     * @return The return value an instance of the target type that has
     * had its value set by stream input from the specified string.
     */
    template <class Type>
    Type
    convertString(std::string const& inputString);


    /**
     * This function returns a single string comprising copies of all of
     * the strings in inputStringVector, interposed by the copies of
     * separator.  For example:
     *
     * @code
     *   joinString(["Hi", "there", "honey", "bear"], "/");
     * @endcode
     *
     * will return "Hi/there/honey/bear".
     *
     * @param inputStringVector This argument specifies the strings to
     * be concatenated.
     *
     * @param separator This argument specifies text to be repeated
     * between consecutive elements of inputStringVector.
     *
     * @return The return value is the concatenation of the elements of
     * inputStringVector, with the contents of separator interposed
     * between elements.
     */
    std::string
    joinString(const std::vector<std::string>& inputStringVector,
               std::string const& separator="");


    /**
     * This function returns a copy of the input argument in which every
     * uppercase character has been replaced with its lowercase
     * equivalent.  For example:
     *
     * @code
     *   makeLowerCase("hI TheRe!")
     * @endcode
     *
     * will return "hi there!"
     *
     * @param inputString This argument specifies the string that will
     * be converted to lower case.
     *
     * @return The return value is the converted string.
     */
    std::string
    makeLowerCase(std::string const& inputString);


    /**
     * This function returns a copy of the input argument in which every
     * lowercase character has been replaced with its uppercase
     * equivalent.  For example:
     *
     * @code
     *   makeUpperCase("hI TheRe!")
     * @endcode
     *
     * will return "HI THERE!"
     *
     * @param inputString This argument specifies the string that will
     * be converted to upper case.
     *
     * @return The return value is the converted string.
     */
    std::string
    makeUpperCase(std::string const& inputString);


    /**
     * This function copies a string, replacing non-overlapping
     * occurrences of a target string with the specified replacement.
     * For example:
     *
     * @code
     *   replaceString("hohoohooohooooho", "oo", "ii")
     * @endcode
     *
     * gives
     *
     *   "hohiihiiohiiiiho"
     *
     * @param inputString This argument species the string that is to
     * be copied.
     *
     * @param target This argument specifies the substring that is to
     * be replaced.
     *
     * @param replacement This argument specifies the string that is to
     * be substituted for each instance of the target.
     *
     * @return The return value is a string in which all non-overlapping
     * instances of target have been replaced with replacement.
     */
    std::string
    replaceString(std::string const& inputString,
                  std::string const& target,
                  std::string const& replacement);


    /**
     * This function divides inputString around instances of delimiter.
     * For example:
     *
     * @code
     *   splitString("04//22/2001", "/")
     * @endcode
     *
     * gives
     *
     *   ["04", "22", "2001"]
     *
     * and
     *
     * @code
     *   splitString("04 - 22 - 2001", " - ")
     * @endcode
     *
     * gives
     *
     *   ["04", "22", "2001"]
     *
     * Setting includeNullStrings to true makes the return value include
     * even strings of zero length.
     *
     * @code
     *   splitString("04//22/2001", "/", true)
     * @endcode
     *
     * gives
     *
     *   ["04", "", "22", "2001"]
     *
     * and
     *
     * @code
     *   splitString("hellohehoundhe", "he", true)
     * @endcode
     *
     * gives
     *
     *   ["", "llo", "hound", ""]
     *
     * @param inputString This argument specifies the string that is to
     * be split.
     *
     * @param delimiter This argument the substring around which
     * inputString should be divided.
     *
     * @param includeNullStrings This argument specifies whether
     * zero-length substrings should be included in the output.
     *
     * @param maxSplit If this argument is non-zero, it limits the
     * number of splits performed.  If maxSplit is non-zero, then the
     * returned vector of strings will have at most (maxSplit + 1)
     * elements, with the first maxSplit elements being formed by
     * splitting chunks off the beginning of the string, and the final
     * element being the remainder of the string after the first
     * maxSplit splits are performed.  If this argument is zero, then
     * the number of splits will not be limited.
     *
     * @return The return value is a vector of substrings as described
     * above.
     */
    std::vector<std::string>
    splitString(std::string const& inputString,
                std::string const& delimiter,
                bool includeNullStrings=false,
                size_t maxSplit=0);


    /**
     * This function removes whitespace from the beginning and end of
     * a string.  If the second argument is specified, it indicates which
     * characters should be considered to be whitespace.
     *
     * @param inputString This is the string to be stripped.
     *
     * @param whiteSpace The characters of this string define what qualifies
     * as whitespace.
     *
     * @return A copy of inputString from which leading and trailing
     * whitespace has been removed.
     */
    std::string
    stripString(std::string const& inputString,
                std::string const& whiteSpace=" \t\n");


    /**
     * This function returns a copy of the input argument in which
     * end-of-line markers have been inserted to wrap the string at a
     * specified line length.  Wrapping will only occur at places where
     * there is whitespace.  This means that long sections with no
     * whitespace will not be broken, and may exceed the requested line
     * length.  For example:
     *
     * @code
     *   wrapString("This is a string to be wrapped.\nYesyesyes indeed it is",
     *              "", 8, " ", "\n")
     * @endcode
     *
     * will return
     *
     *   "This is\na string\nto be\nwrapped.\nYesyesyes\nindeed\nit is"
     *
     * and
     *
     * @code
     *   wrapString("This is a string to be wrapped.\nYesyesyes indeed",
     *              "> ", 10, " ", "\n")
     * @endcode
     *
     * will return
     *
     *   "This is a\n> string\n> to be\n> wrapped.\n> Yesyesyes\n> indeed"
     *
     * @param inputString This argument is the string to be wrapped.
     *
     * @param fillPrefix This prefix will be inserted immediately after
     * each newline in the output string.
     *
     * @param width This argument specifies the target line length for
     * the returned string.  Individual lines will only exceed this
     * length if there are words in the input string that have more
     * than width characters.
     *
     * @param whitespace This argument specifies which characters should
     * count as whitespace when breaking inputString.  The string may
     * be broken around an character present in whitespace.
     *
     * @param eolString This argument specifies the end-of-line marker.
     * This is used for two things: its presence is used to identify
     * existing newlines in the input string, and it is inserted into
     * the output string to introduce a newline.
     *
     * @return The return value is a wrapped string.
     */
    std::string
    wrapString(std::string const& inputString,
               std::string const& fillPrefix="",
               size_t width=78,
               std::string const& whitespace=" \t\n",
               std::string const& eolString="\n");

  } // namespace utilities

} // namespace brick


/* ******************************************************************
 * Definitions of inline and template functions.
 * ***************************************************************** */

namespace brick {

  namespace utilities {

    // Converts the string inputString to Type.
    template <class Type>
    Type
    convertString(std::string const& inputString)
    {
      Type returnValue;
      std::istringstream inputStream(inputString);
      inputStream >> returnValue;
      if(!inputStream) {
        std::ostringstream message;
        message << "Can't convert malformed input string: " << inputString
                << std::endl;
        BRICK_THROW(ConversionException, "convertString()",
                    message.str().c_str());
      }
      return returnValue;
    }

  } // namespace utilities

} // namespace brick

#endif /* #ifndef BRICK_UTILITIES_STRINGMANIPULATION_HH */
