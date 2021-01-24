#!/usr/bin/python

import os
import string
import sys

stateIndex = 0

STATE_LOOKING_FOR_OPENING_COMMENT = stateIndex
stateIndex += 1

STATE_LOOKING_FOR_EMAIL_ADDRESS = stateIndex
stateIndex += 1

STATE_LOOKING_FOR_COMMENT_END = stateIndex
stateIndex += 1

STATE_FINISHED = stateIndex
stateIndex += 1


oldEmailAddress = 'dlr@cs.cmu.edu'
newEmailAddress = 'dlr@davidlarose.com'

if __name__ == '__main__':

  # Read the input
  inputFileName = sys.argv[1]
  inputFile = open(inputFileName)
  inputLines = inputFile.readlines()
  inputFile.close()

  # If we made it this far, it's time to do a backup.
  backupFileName = inputFileName + '.ffhBackup'
  if os.path.isfile(backupFileName):
    raise IOError('Backup file %s already exists.' % backupFileName)
  # end if
  # os.copyfile(inputFileName, backupFileName)
  # os.system('cp -a %s %s' % (inputFileName, backupFileName))
  os.rename(inputFileName, backupFileName)

  # Overwrite the input file.
  outputFile = open(inputFileName, 'w')
  
  # Now start parsing file
  currentState = STATE_LOOKING_FOR_OPENING_COMMENT
  
  for line in inputLines:

    if currentState == STATE_LOOKING_FOR_OPENING_COMMENT:
      if line.strip()[:2] == '/*':
        currentState = STATE_LOOKING_FOR_EMAIL_ADDRESS
      # end if
      outputFile.write(line)
      continue
    # end if

    if currentState == STATE_LOOKING_FOR_EMAIL_ADDRESS:
      if line.find(oldEmailAddress) > 0:
        currentState = STATE_FINISHED
        line = line.replace(oldEmailAddress, newEmailAddress)
      # end if
      outputFile.write(line)
      continue
    # end if

    if currentState == STATE_FINISHED:
      outputFile.write(line)
      continue
    # end if

  # end for

  outputFile.close()

  if currentState != STATE_FINISHED:
    raise IOError('Trouble processing input file %s.' % inputFileName)
  # end if
# end if
