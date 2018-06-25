/**
***************************************************************************
* @file stringManipulation.cc
*
* Source file defining routines for working with strings.
*
* Copyright (C) 2003-2011, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

// #include <ctype.h>
#include <algorithm>
#include <string>
#include <brick/utilities/stringManipulation.hh>

namespace brick {

  namespace utilities {

    // This function takes an input string and returns a string in which
    // all special shell characters have been escaped.
    std::string
    cleanString(std::string const& inputString,
                std::string const& specialCharacters,
                char quoteCharacter,
                bool alreadyQuoted)
    {
      std::string outputString;
      std::string::size_type index = 0;
      while(index < inputString.size()) {
        // Look for the next special character.
        std::string::size_type nextIndex =
          inputString.find_first_of(specialCharacters, index);

        // Add all of the normal characters before the next special
        // character to the output string.
        if(nextIndex != index) {
          outputString +=
            inputString.substr(index, nextIndex - index);
        }

        // If we didn't find a special character, we're done.
        if(nextIndex == std::string::npos) {
          break;
        }

        // If we did find a special character, quote it unless it's a
        // quote character itself and we've been told not to quote quote
        // characters.
        if(alreadyQuoted && (inputString[nextIndex] == quoteCharacter)) {
          // It's a quote character, and we're not supposed to quote it,
          // so pass it on unchanged.
          outputString += inputString[nextIndex];

          // And pass the character that it's quoting, too.
          ++nextIndex;
          if(nextIndex < inputString.size()) {
            outputString += inputString[nextIndex];
          }
        } else {
          // We're allowed to quote this character.
          outputString += quoteCharacter;
          outputString += inputString[nextIndex];
        }

        // And start the next search from one past the last character we
        // added.
        index = nextIndex + 1;
      }

      // Return the escaped version of the string.
      return outputString;
    }


    // This function returns a single string comprising copies of all of
    // the strings in inputStringVector, interposed by the copies of
    // separator.
    std::string
    joinString(const std::vector<std::string>& inputStringVector,
               std::string const& separator)
    {
      std::ostringstream outputStream;
      std::vector<std::string>::const_iterator inputIterator =
        inputStringVector.begin();

      if(inputIterator != inputStringVector.end()) {
        outputStream << *inputIterator;
        ++inputIterator;
      }

      while(inputIterator != inputStringVector.end()) {
        outputStream << separator << *inputIterator;
        ++inputIterator;
      }
      return outputStream.str();
    }


    // This function returns a copy of the input argument in which every
    // uppercase character has been replaced with its lowercase
    // equivalent.
    std::string
    makeLowerCase(std::string const& inputString)
    {
      std::string result(inputString.size(), '\0');
      std::transform(inputString.begin(), inputString.end(),
                     result.begin(), ::tolower);
      return result;
    }


    // This function returns a copy of the input argument in which every
    // lowercase character has been replaced with its uppercase
    // equivalent.
    std::string
    makeUpperCase(std::string const& inputString)
    {
      std::string result(inputString.size(), '\0');
      std::transform(inputString.begin(), inputString.end(),
                     result.begin(), ::toupper);
      return result;
    }


    // Copies inputString, replacing non-overlapping occurrences of
    // target with replacement.
    std::string
    replaceString(std::string const& inputString,
                  std::string const& target,
                  std::string const& replacement)
    {
      std::string outputString;
      std::string::size_type startIndex = 0;
      while(startIndex <= inputString.size()) {
        std::string::size_type targetPosition =
          inputString.find(target, startIndex);
        if(targetPosition == std::string::npos) {
          // target not found, so add the rest of inputString & quit.
          outputString += inputString.substr(startIndex);
          break;
        }
        // Target found, so add up to it's beginning.
        outputString +=
          inputString.substr(startIndex, targetPosition - startIndex);
        // And then add replacement.
        outputString += replacement;
        // Get ready to start searching again.
        startIndex = targetPosition + target.size();
      }
      return outputString;
    }


    // Divides inputString around instances of delimiter.
    std::vector<std::string>
    splitString(std::string const& inputString,
                std::string const& delimiter,
                bool includeNullStrings,
                size_t maxSplit)
    {
      std::vector<std::string> stringParts;
      std::string::size_type startIndex = 0;
      size_t splitCount = 0;
      while(startIndex < inputString.size()) {
        std::string::size_type endIndex =
          inputString.find(delimiter, startIndex);
        if(endIndex == std::string::npos) {
          // npos is hard to work with.  This is more convenient.
          endIndex = inputString.size();
        }
        // Terminate early if we've exceeded a non-zero maxSplit.
        if(maxSplit != 0 && splitCount >= maxSplit) {
          endIndex = inputString.size();
        }
        // Delimiter not at very beginning of remaining string?
        if(includeNullStrings || (endIndex != startIndex)) {
          // Push the next section
          stringParts.push_back(inputString.substr(startIndex,
                                                   endIndex - startIndex));
          ++splitCount;
        }
        startIndex = endIndex + delimiter.size();
      }
      if(includeNullStrings && (startIndex == inputString.size())) {
        stringParts.push_back(std::string(""));
      }
      return stringParts;
    }


    // This function removes whitespace from the beginning and end of
    // a string.
    std::string
    stripString(std::string const& inputString, std::string const& whiteSpace)
    {
      std::string::size_type startIndex =
        inputString.find_first_not_of(whiteSpace);
      std::string::size_type endIndex =
        inputString.find_last_not_of(whiteSpace);
      if((startIndex == std::string::npos) || (endIndex == std::string::npos)) {
        return std::string();
      }
      return inputString.substr(startIndex, endIndex - startIndex + 1);
    }


    // This function returns a copy of the input argument in which
    // end-of-line markers have been inserted to wrap the string at a
    // specified line length.
    std::string
    wrapString(std::string const& inputString,
               std::string const& fillPrefix,
               size_t width,
               std::string const& whitespace,
               std::string const& eolString)
    {
      enum WrapState {
        BRICK_WS_STARTING_LINE,
        BRICK_WS_SEARCHING_FOR_EOL,
        BRICK_WS_SEARCHING_FOR_WHITESPACE0,
        BRICK_WS_SEARCHING_FOR_WHITESPACE1,
        BRICK_WS_EATING_BEGINNING_WHITESPACE,
        BRICK_WS_EATING_ENDING_WHITESPACE,
        BRICK_WS_FINISHED
      };

      // We'll use an ostringstream to format the result.
      std::ostringstream outputStream;

      // The position in inputString at which we'll start searching for
      // whitespace & eol characters.
      std::string::size_type startIndex = 0;

      // The first column of the current line that is available for
      // characters from inputString.  This is always zero for the
      // first line, and fillPrefix.size() for subsequent lines.
      std::string::size_type startColumn = 0;

      // Whenever this variable is greater than startIndex, it will
      // indicate the position of the next newline in the inputString.
      std::string::size_type eolIndex = 0;

      // Whenever this variable is greater than startIndex, it will
      // indicate the position of the next whitespace in the inputString.
      std::string::size_type whitespaceIndex0 = 0;

      // Whenever this variable is greater than whitespaceIndex0, it
      // will indicate the position of the next subsequent whitespace in
      // the inputString.
      std::string::size_type whitespaceIndex1 = 0;

      // Whenever this variable is greater than startIndex, it will
      // indicate the position of the next non-whitespace character.
      std::string::size_type nonWhitespaceIndex = 0;

      // Loop until done processing the entire string.
      WrapState currentState = BRICK_WS_STARTING_LINE;
      while(currentState != BRICK_WS_FINISHED) {

        // Sanity check here to make sure we don't run off the end of
        // the string.  This simplifies the code below a little.
        if(startIndex >= inputString.size()) {
          currentState = BRICK_WS_FINISHED;
          break;
        }

        // This switch implements the parsing algorithm.
        switch(currentState) {
        case BRICK_WS_STARTING_LINE:
          if(inputString.size() - startIndex <= width - startColumn) {
            // Remaining text will all fit on one line.  We're done!
            outputStream << inputString.substr(startIndex, std::string::npos);
            currentState = BRICK_WS_FINISHED;
          } else {
            // String doesn't fit.  Look for line breaks.
            currentState = BRICK_WS_SEARCHING_FOR_EOL;
          }
          break;
        case BRICK_WS_SEARCHING_FOR_EOL:
          // Look for existing line breaks.
          if(eolIndex <= startIndex) {
            eolIndex = inputString.find(eolString, startIndex);
          }
          if(eolIndex - startIndex <= width - startColumn) {
            // There's a line break early enough in the string that we
            // can just break there without exceeding width.
            outputStream
              << inputString.substr(startIndex, eolIndex - startIndex)
              << eolString;
            startIndex = eolIndex + eolString.size();

            // We just printed a line break.  Eat whitespace so we don't
            // start a line with spaces, and so we're sure there remains
            // non-whitespace to print.
            currentState = BRICK_WS_EATING_BEGINNING_WHITESPACE;
          } else {
            // No appropriate line break found.  Look for whitespace at
            // which to break.
            currentState = BRICK_WS_SEARCHING_FOR_WHITESPACE0;
          }
          break;
        case BRICK_WS_SEARCHING_FOR_WHITESPACE0:
          // Search for next whitespace.
          if(whitespaceIndex0 <= startIndex) {
            whitespaceIndex0 = inputString.find_first_of(whitespace, startIndex);
          }
          if(whitespaceIndex0 == std::string::npos) {
            // No whitespace available.  Nothing we can do.  Simply
            // write out the rest of the string.
            outputStream << inputString.substr(startIndex, std::string::npos);
            currentState = BRICK_WS_FINISHED;
          } else if(whitespaceIndex0 - startIndex > width - startColumn) {
            // Next whitespace is far enough away that we exceed width.
            // We're stuck with an overlong line.
            outputStream
              << inputString.substr(startIndex, whitespaceIndex0 - startIndex);
            startIndex = whitespaceIndex0 + 1;

            // Before we print a linebreak, we need to make sure that
            // there's more non-whitespace to print.  In other words, we
            // need to be sure that we actually need another line of
            // output.
            currentState = BRICK_WS_EATING_ENDING_WHITESPACE;
          } else {
            // Found whitespace that's early enough to break without
            // exceeding width, but is there more whitespace that's even
            // better?
            currentState = BRICK_WS_SEARCHING_FOR_WHITESPACE1;
          }
          break;
        case BRICK_WS_SEARCHING_FOR_WHITESPACE1:
          // Make sure we won't run off the end of the string by
          // searching for more whitespace.
          if(whitespaceIndex0 == inputString.size() - 1) {
            // Looks like we've reached the end of the string.  Write it out.
            outputStream
              << inputString.substr(startIndex, whitespaceIndex0 - startIndex);
            currentState = BRICK_WS_FINISHED;
          } else {
            // Search for next whitespace without forgetting the current
            // whitespace.
            whitespaceIndex1 = whitespaceIndex0;
            whitespaceIndex0 =
              inputString.find_first_of(whitespace, whitespaceIndex1 + 1);
            if(whitespaceIndex0 - startIndex <= width - startColumn) {
              // Found whitespace that's early enough to break at without
              // exceeding width, but is there another whitespace that's
              // even better?
              currentState = BRICK_WS_SEARCHING_FOR_WHITESPACE1;
            } else {
              // No more whitespace before we exceed width.  Break the
              // string at the best position found.
              outputStream
                << inputString.substr(startIndex, whitespaceIndex1 - startIndex);
              startIndex = whitespaceIndex1 + 1;

              // Before we print a linebreak, we need to make sure that
              // there's more non-whitespace to print.  In other words, we
              // need to be sure that we actually need another line of
              // output.
              currentState = BRICK_WS_EATING_ENDING_WHITESPACE;
            }
          }
          break;
        case BRICK_WS_EATING_BEGINNING_WHITESPACE:
          // This state happens immediately after we copy a line break
          // from the input string to the output.
          //
          // Find next non-whitespace character.
          if(nonWhitespaceIndex <= startIndex) {
            nonWhitespaceIndex =
              inputString.find_first_not_of(whitespace, startIndex);
          }
          if(nonWhitespaceIndex == std::string::npos) {
            // No more non-whitespace.  We're done!
            currentState = BRICK_WS_FINISHED;
          } else {
            // Found non-whitespace.  Prepare for next line.
            outputStream << fillPrefix;
            startIndex = nonWhitespaceIndex;
            startColumn = fillPrefix.size();
            currentState = BRICK_WS_STARTING_LINE;
          }
          break;
        case BRICK_WS_EATING_ENDING_WHITESPACE:
          // This state happens immediately before we print a line break
          // that wasn't part of the original string.
          //
          // Find next non-whitespace character.
          if(nonWhitespaceIndex <= startIndex) {
            nonWhitespaceIndex =
              inputString.find_first_not_of(whitespace, startIndex);
          }
          if(nonWhitespaceIndex == std::string::npos) {
            // No more non-whitespace.  We're done!
            currentState = BRICK_WS_FINISHED;
          } else {
            // Found non-whitespace.  We'll need another line.
            outputStream << eolString << fillPrefix;
            startIndex = nonWhitespaceIndex;
            startColumn = fillPrefix.size();
            currentState = BRICK_WS_STARTING_LINE;
          }
          break;
        default:
          // Will never get here.
          break;
        }
      }

      // Done.  Return the string we just formatted.
      return outputStream.str();
    }

  } // namespace utilities

} // namespace brick
