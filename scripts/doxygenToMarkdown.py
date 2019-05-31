#!/usr/bin/python

import re
import sys

def markdownify(inputFileName, outputFileName):
  # Read the input file
  inputFile = open(inputFileName)
  inputData = inputFile.read()
  inputFile.close()

  # Discard initial comment characters.
  inputData = re.sub(r'^/\*\*', '', inputData)

  # Discard doxygen file identifier.
  inputData = re.sub(r'@file .*\n', '', inputData)

  # Change doxygen title to markdown header1.
  inputData = re.sub(r'@mainpage', '#', inputData)

  # Change doxygen section to markdown header2.
  inputData = re.sub(r'@section( \w+_sec)?', '##', inputData)

  # Change doxygen subsection to markdown header3.
  inputData = re.sub(r'@subsection', '###', inputData)

  # Change doxygen subsubsection to markdown header4.
  inputData = re.sub(r'@subsubsection', '####', inputData)

  # Change doxygen verbatim sections to markdown code.
  inputData = re.sub(r'@(end)?verbatim', '```', inputData)

  # Change doxygen code sections to markdown code.
  inputData = re.sub(r'@(end)?code', '```', inputData)

  # Remove closing comment characters.
  inputData = re.sub(r'\*\*/\n?', '', inputData)

  # Save resulting markdown file.
  outputFile = open(outputFileName, 'w')
  outputFile.write(inputData)
  outputFile.close()
# end def

if __name__ == '__main__':
  markdownify(sys.argv[1], sys.argv[2])
# end if

