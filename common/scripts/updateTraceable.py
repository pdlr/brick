#!/usr/bin/python

## This file is intended to simplify the process of adding
## BEGIN_TRACEABLE/END_TRACEABLE macro pairs to C++ source code.  To
## use it, start by editing your C++ file so that every function you'd
## like to instrument has a BEGIN_TRACEABLE macro.  Don't worry about
## putting in the END_TRACEABLE macros.  Then, run this file with the
## filename of your C++ file as the only argument:
##
##   python updateTraceable.py mySourceFile.cpp
##
## Then go back and check the file to make sure the END_TRACEABLE
## macros have been added to your satisfaction.

import os
import re
import shutil
import string
import sys

## This function takes details about the current function and formats
## the corresponing END_TRACEABLE macro text.  For example,
##
##   allArgumentsFuncor('foo', '(int bar, int baz)', '  ')
##
## will return
##
##   '  END_TRACEABLE(foo, (bar, baz));'
##
def allArgumentsFunctor(functionName, argumentList, indent):
  ## We'll break lines so they don't wrap on 80 column screens.
  maxLineLength = 79

  ## Parse the argument list to extract the names of the arguments.
  argumentNameMatches = re.finditer(r'(\w*)\s*[,\)]', argumentList)
  argumentNameList = map(lambda x: x.group(1), argumentNameMatches)

  ## Start formatting the macro.
  initialText = indent + 'END_TRACEABLE('
  macroTextComponents = [initialText]

  ## Now we have a decision to make.  If the function name isn't too long,
  ## we'll stay on the same line.  If it is too long to fit on this line,
  ## we'll drop down to the next line and increase the indent.
  if (len(initialText) + len(functionName) + 1) >= maxLineLength:
    localIndent0 = indent + '  '
    macroTextComponents.append('\n')
    macroTextComponents.append(localIndent0)
  else:
    localIndent0 = (
      indent + string.join([' '] * (len(initialText) - len(indent)), ''))
  # end if

  ## Now add the function name.
  macroTextComponents.append(functionName)
  macroTextComponents.append(', ')

  ## Figure out how many columns we'll need for the argument list.
  argumentsLength = len(string.join(argumentNameList, ', ')) + 2

  ## Check whether the remaining text will fit on the same line.
  if ((len(localIndent0) + len(functionName) + argumentsLength + 1)
      < maxLineLength):
    ## Yes!  Finish up the macro text.
    macroTextComponents.append('(')
    macroTextComponents.append(string.join(argumentNameList, ', '))
    macroTextComponents.append('));')
  else:
    ## No!  Insert carriage returns and indent as needed.
    macroTextComponents.append('\n')
    macroTextComponents.append(localIndent0)
    macroTextComponents.append('(')

    ## Insert the first argument, if necessary.
    startingPosition = len(localIndent0) + 1
    currentPosition = startingPosition
    if len(argumentNameList) >= 1:
      macroTextComponents.append(argumentNameList[0])
      currentPosition = startingPosition + len(argumentNameList[0])
    # end if

    ## Continue, inserting subsequent arguments.
    if len(argumentNameList) >= 2:
      for argumentName in argumentNameList[1:]:
        ## Before inserting each argument, check to see if a carriage
        ## return is needed.
        if (currentPosition + len(argumentName) + 4) >= maxLineLength:
          macroTextComponents.append(',\n' + localIndent0)
          currentPosition = startingPosition
        else:
          macroTextComponents.append(', ')
          currentPosition = currentPosition + 2
        # end if
        macroTextComponents.append(argumentName)
        currentPosition = currentPosition + len(argumentName)
      # end for
    # end if
    macroTextComponents.append('));')
  # end if

  return string.join(macroTextComponents, '')
# end def


## This class does all the work of adding macros.
class TraceableParser:
  def __init__(self):

    ## The regular expression for matching the BEGIN_TRACEABLE macro
    ## should be run with MULTILINE turned on, and breaks down as follows:
    ##   ^        - Matches the beginning of the line.
    ##   ([ \t]*) - Matches the indent preceding the macro and saves it.
    ##   BEGIN_TRACEABLE - Matches the macro itself.
    ##   ;?       - Matches any trailing semicolon.
    ##   \s*      - Matches any whitespace before the end of the line.
    ##   \n?      - Matches the trailing '\n', if present.
    self.m_beginMacroPattern = r'^([ \t]*)BEGIN_TRACEABLE;?\s*?\n?'

    ## The regular expression for matching the BEGIN_TRACEABLE macro
    ## should be run with MULTILINE and DOTALL turned on, and breaks
    ## down as follows:
    ##   ^      - Matches the beginning of the line.
    ##   [ \t]* - Matches the indent preceding the macro.
    ##   END_TRACEABLE - Matches the macro itself.
    ##   \(.*?\)+;? - Matches the argument list, if present.
    ##   ;?     - Matches any trailing semicolon.
    ##   \s*    - Matches any whitespace before the end of the line.
    ##   \n?    - Matches the trailing '\n', if present.
    self.m_endMacroPattern = r'^[ \t]*END_TRACEABLE\(.*?\)+;?\s*?\n?'
    
    ## The regular expression for matching function prototypes should
    ## be run with MULTILINE and DOTALL turned on, and breaks down as
    ## follows:
    ##   [~,:\s] - Matches only if the first characer is not ',', ':', or
    ##             whitespace.  This eliminates matches with constructor
    ##             initialization list entries.  Note that this means we
    ##             won't match a function definition if its first character
    ##             is the first non-whitespace character in the file.  In
    ##             our code this never happens.
    ##   \s+     - Matches any whitespace preceding the function/class name.
    ##             Note that we require there to be some whitespace.  This
    ##             is a departure from the C++ grammar, but it makes the
    ##             regex easier to write.
    ##   ((?:\w+\s*::\s*)?\w+) - Matches the function name and saves the match
    ##                           as a group.  The matching is as follows:
    ##     (?:\w+\s*::\s*)? - Matches class name qualifiers, like "MyClass::".
    ##                        The trailing '?' makes it optional, and the
    ##                        "?:" after the opening parenthesis prevents
    ##                        the match from counting as a group.  The contents
    ##                        match as follows:
    ##       \w+  - Matches the class name.
    ##       \s*  - Matches any whitespaces between the class name and the
    ##              "::".
    ##       ::   - Matches "::"
    ##       \s*  - Matches any whitespaces between the "::" and the function
    ##              name.
    ##     \s*           - Matches any whitespace between the class name
    ##                     and the function name.
    ##     ~?        - Matches the '~' in destructor names. 
    ##     \w+       - Matches the function name, like "foo".
    ##   \s*       - Matches any whitespace between the function name and
    ##               argument list.
    ##   (\([^{}]*?\)) - Matches the argument list and saves the match as a
    ##                   group.
    ##
    ## We specify the DOTALL flag so that '.' will match newline characters
    ## as well as all others.
    self.m_prototypePattern = (
      r'[^,:\s]\s+((?:\w+\s*::\s*)?~?\w+)\s*(\([^{}]*?\))')

    ## The macro for matching a function prototype + BEGIN_TRACEABLE
    ## macro is simply the conjunction of the two previously defined
    ## regular expressions, allowing for free text in between.
    self.m_prototypeAndBeginMacroPattern = (
      self.m_prototypePattern + r'[^}]*?' + self.m_beginMacroPattern)


    ## Compile the various patterns.
    reFlags = re.MULTILINE | re.DOTALL
    self.m_beginMacroRegEx = re.compile(self.m_beginMacroPattern, reFlags)
    self.m_endMacroRegEx = re.compile(self.m_endMacroPattern, reFlags)
    self.m_prototypeRegEx = re.compile(self.m_prototypePattern, reFlags)
    self.m_prototypeAndBeginMacroRegEx = re.compile(
      self.m_prototypeAndBeginMacroPattern, reFlags)
  # end def


  ## This method searches for BEGIN_TRACEABLE macros, and inserts
  ## corresponding END_TRACEABLE macros.
  def addEndMacros(self, endMacroFunctor):
    # Find all instances of the begin macro text and associated function
    # prototypes.
    beginMacroMatchIter = self.m_prototypeAndBeginMacroRegEx.finditer(
      self.m_fileText)

    # We'll process the file from end to start so that as we
    # process each line we don't invalidate the following ones.
    beginMacroMatches = []
    for macroMatch in beginMacroMatchIter:
      beginMacroMatches.append(macroMatch)
    # end for
    beginMacroMatches.reverse()

    ## Consider each instance of BEGIN_TRACEABLE in turn.
    for macroMatch in beginMacroMatches:
      ## Figure out where the corresponding END_TRACEABLE should go.
      scopeEnd = self.findScopeEnd(macroMatch.end(), '{', '}', 1)
      insertionPoint = self.rewindOverIndent(scopeEnd)

      ## Recover information that should go into the END_TRACEABLE macro
      ## arguments.
      functionName = macroMatch.group(1)
      argumentList = macroMatch.group(2)
      prettyFunctionName = re.sub(r'\s+', '', functionName)

      ## Format the END_TRACEABLE macro using the supplied function.
      macroText = endMacroFunctor(
        prettyFunctionName, argumentList, macroMatch.group(3))

      ## And insert the text.
      self.insertText(macroText + '\n', insertionPoint)
    # end for
  # end def


  ## This method removes all existing END_TRACEABLE macros, to make
  ## way for auto-generated ones.
  def deleteEndMacros(self):
    self.m_fileText = self.m_endMacroRegEx.sub("", self.m_fileText)
  # end def


  ## This method loads a file in preparation for parsing.
  def load(self, fileName):
    inputFile = open(fileName)
    self.m_fileText = string.join(inputFile.readlines(), '')
    inputFile.close()
  # end def


  ## This method saves the processed text to a file.
  def save(self, fileName):
    outputFile = open(fileName, 'w')
    outputFile.write(self.m_fileText)
    outputFile.close()
  # end def


  # ============ Private methods follow ============


  ## This method searches from the specified position (which is an
  ## index into m_fileText) to find the end of the enclosing scope.
  ## The arguments openingCharacter and closingCharacter define the
  ## scope delimiters.  Argument scopeLevel is used to set the
  ## starting scope level.  This is useful if, for example, the search
  ## starts inside and if clause, but you want to keep searching and
  ## find the end of the enclosing scope.  In this case you'd set
  ## scopeLevel to 1, indicating that there is one level of scope
  ## level to be ignored.  
  def findScopeEnd(self, position, openingCharacter, closingCharacter,
                   scopeLevel):
    ## Iterate through the file looking for opening and closing
    ## characters until we get down to scopeLevel == 0.
    closingPosition = position;
    while closingPosition < len(self.m_fileText):
      if self.m_fileText[closingPosition] == openingCharacter:
        scopeLevel = scopeLevel + 1
      # end if
      if self.m_fileText[closingPosition] == closingCharacter:
        scopeLevel = scopeLevel - 1
      # end if
      if scopeLevel <= 0:
        return closingPosition
      # end if
      closingPosition = closingPosition + 1
    # end while
    raise ValueException, ('Can\'t find closing \'%s\' at position %d.'
                           % (closingCharacter, position))
  # end def


  ## This method simply inserts the specified text at the specified
  ## point in m_fileText.
  def insertText(self, macroText, insertionPoint):
    self.m_fileText = (
      self.m_fileText[:insertionPoint] + macroText
      + self.m_fileText[insertionPoint:])
  # end def


  ## This method is used to find the beginning of the line at which
  ## the END_TRACEABLE macro should be inserted.
  def rewindOverIndent(self, position):
    ## Begin by handling the trivial case.
    if position <= 0:
      return 0;
    # end if

    ## Back up until we see a character that's not a space or a tab.
    position = position - 1
    while position >= 0:
      if self.m_fileText[position] not in [' ', '\t']:
        position = position + 1
        break
      # end if
      position = position - 1
    # end while
    return position
  # end def
  
# end class


## This function makes a copy of the original file, in case of errors.
def backupInputFile(inputFileName):
  backupFileName = inputFileName + '.utBackup'
  if os.path.isfile(backupFileName):
    raise IOError, ('Backup file %s already exists.' % backupFileName)
  # end if
  shutil.copyfile(inputFileName, backupFileName)
# end def


## This is the usage function.
def usage(argv):
  print 'Usage: %s sourceFile' % argv[0]

  print '  This function parses a C++ source file, and inserts\n'
  print '  END_TRACEABLE macros in any function which contains a\n'
  print '  BEGIN_TRACEABLE macro.'
# end def


## Main routine.
if __name__ == '__main__':
  if len(sys.argv) != 2:
    usage(sys.argv)
    sys.exit(0)
  # end if

  backupInputFile(sys.argv[1])
  parser = TraceableParser()
  parser.load(sys.argv[1])
  parser.deleteEndMacros()
  parser.addEndMacros(allArgumentsFunctor)
  parser.save(sys.argv[1])
# end if
