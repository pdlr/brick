#!/usr/bin/python

import os
import string
import sys


if __name__ == '__main__':

  # Make sure we can read the input
  inputFileName = sys.argv[1]
  inputFile = open(inputFileName)
  # inputLines = inputFile.readlines()
  inputFile.close()

  # Make sure we can open the output file.
  outputFileName = sys.argv[2]
  outputFile = open(outputFileName, 'w')
  outputFile.close()

  # We'll need a temporary file to work with.
  tempFileName = outputFileName + '.deleteMe'
  
  # Copy the data from inputFile once and for all.
  os.system("cp %s %s" % (inputFileName, outputFileName))

  # Discard CVS-style automatically substituted keywords.
  os.system("sed "
            + r"-e '/^.*\$Revision:.*\$/ d' "
            + r"-e '/^.*\$Date:.*\$/ d' "
            + "< %s > %s" % (outputFileName, tempFileName))
  os.system("cp %s %s" % (tempFileName, outputFileName))

  # Fix include guards and macro names.
  #  - 1st and 2nd commands remove leading and trailing "_" from DLR
  #    macros.
  #  - 3rd command turns macros that start with "DLR" (no underscore)
  #    into macros that start with "DLR_".
  #  - 4th command turns macros that start with "DLR_" (no underscore)
  #    into macros that start with "BRICK_".
  #  - 5th command turns macros that start with "_H" into macros that
  #    start with "_HH".
  os.system("sed "
            + r"-e 's/_DLR/DLR/' "
            + r"-e 's/\(DLR[^ 	]*\)_$/\1/' "
            + r"-e 's/\(DLR[^ 	]*\)_\([ 	]\)/\1\2/' "
            + r"-e 's/DLR\([^_ 	]\)/DLR_\1/' "
            + r"-e 's/DLR_/BRICK_/' "
            + r"-e 's/_H/_HH/' "
            + "< %s > %s" % (outputFileName, tempFileName))
  os.system("cp %s %s" % (tempFileName, outputFileName))
  
  # Fix include directives and doxygen file comments.
  os.system("sed "
            + "-e 's/dlrCommon\//brick\/common\//' "
            + "-e 's/\.h$/\.hh/' "
            + "-e 's/\.h\>/\.hh\>/' "
            + "< %s > %s" % (outputFileName, tempFileName))
  os.system("cp %s %s" % (tempFileName, outputFileName))

  # Fix namespace declarations.
  os.system("sed "
            + "-e 's/namespace dlr/namespace brick/' "
            + "< %s > %s" % (outputFileName, tempFileName))
  os.system("cp %s %s" % (tempFileName, outputFileName))

  # Fix namespace qualifiers.
  os.system("sed "
            + "-e 's/dlr::/brick::/' "
            + "< %s > %s" % (outputFileName, tempFileName))
  os.system("cp %s %s" % (tempFileName, outputFileName))
# end if

