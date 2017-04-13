#!/usr/bin/python

import glob
import os
import re
import shutil
import string
import subprocess
import sys
import tempfile

CMAKE = '/usr/bin/cmake'
GREP = '/bin/grep'
HG = '/usr/bin/hg'
LS = '/bin/ls'
MAKE = '/usr/bin/make'
TAR = '/bin/tar'
# ZIP = '/usr/bin/zip'

# Parse 'hg status' output to see if there are changes to the repo.
MODIFICATION_RE = re.compile('^M|\? .', re.MULTILINE)


def announce(message):
  print message
# end def


def runCommand(args):
  print 'Running command: %s' % string.join(args)
  return subprocess.check_output(args)
# end def


def buildPackage(sourceDir, version, postmortem=False):
  announce('=== Building package ===')

  # On with the show.
  os.chdir(sourceDir)

  print 'Testing for build directories.'
  buildDirList = glob.glob('build*')
  if len(buildDirList) != 0:
    raise IOError('Please delete all build files.')
  # end if
  
  try:
    runCommand([GREP, version, 'VERSION.TXT'])
  except:
    raise IOError('Version number doesn\'t match VERSION.TXT.')
  # end if

  statusString = runCommand([HG, 'status'])
  if MODIFICATION_RE.search(statusString) is not None:
    raise IOError('Repo is not clean.')
  # end if

  buildDir = tempfile.mkdtemp(prefix='brick_build_')
  announce('Build directory is %s' % buildDir)

  try:
    os.chdir(buildDir)

    runCommand([CMAKE, '-D', ('CPACK_PACKAGE_VERSION=%s' % version), sourceDir])
    runCommand([MAKE])
    runCommand([MAKE, 'test'])
    runCommand([MAKE, 'package_source'])

    packageNameBase = 'Brick-%s-Source' % version
    packageNames = [packageNameBase + '.tar.gz',
                    packageNameBase + '.tar.bz2',
                    packageNameBase + '.tar.xz']
    for packageName in packageNames:
      shutil.copyfile(packageName, os.path.join(sourceDir, packageName))
      # end for

    os.chdir(sourceDir)
    shutil.rmtree(buildDir)
  except:
    # We may want to leave the build dir for debugging.
    if not postmortem:
      shutil.rmtree(buildDir)
    # end if
    raise
  # end try

  return packageNames[0]
# end def


def sanitizeVersionInput(version):
  # Get ready for things to go wrong.  This is to avoid repeating code
  # below.
  complaintString = 'Invalid version number: %s' % version

  # We're expecting version to look like '1.11.4' or similar.  Split
  # into major, minor, and patch number.
  components = version.split('.')
  if len(components) != 3:
    raise ValueError(complaintString)
  # end if

  major_minor_patch = (
    int(components[0]), int(components[1]), int(components[2]))

  try:
    sanitizedVersion = '%d.%d.%d' % major_minor_patch
  except:
    raise ValueError(complaintString)
  # end try

  return sanitizedVersion
# end def


def tagRepository(sourceDir, version):
  announce('=== Tagging repository ===')
  
  # We need a civilized version of the revision number to use as a tag
  # in the VCS.
  revisionTag = 'brickLibs_version_%s' % string.replace(version, '.', '_')

  os.chdir(sourceDir)
  runCommand([HG, 'tag', '-f', revisionTag])
# end def


def testPackage(packageName, postmortem=True):
  announce('=== Testing %s ===' % packageName)
  
  testDir = tempfile.mkdtemp(prefix='brick_test_')
  announce('Test directory is %s' % testDir)
  
  try:
    os.chdir(testDir)
    runCommand([TAR, '-zxf', packageName])
    
    buildDir = os.path.join(testDir, 'build')
    os.mkdir(buildDir)
    os.chdir(buildDir)

    fileName = os.path.split(packageName)[1]    # Brick-2.0.0-Source.tar.gz
    baseName = os.path.splitext(fileName)[0]    # Brick-2.0.0-Source.tar
    baseName = os.path.splitext(baseName)[0]    # Brick-2.0.0-Source
    sourceDir = os.path.join(testDir, baseName)
    runCommand([CMAKE, sourceDir])
    runCommand([MAKE])
    runCommand([MAKE, 'test'])
    shutil.rmtree(testDir)
  except:
    # We may want to leave the build dir for debugging.
    if not postmortem:
      shutil.rmtree(testDir)
    # end if
    raise
  # end try
# end def


if __name__ == '__main__':
  if len(sys.argv) != 2:
    print 'Usage: %s version' % sys.argv[0]
    print 'Example: %s 1.11.2' % sys.argv[0]
    sys.exit(65)
  # end if

  version = sys.argv[1]
  sourceDir = os.path.realpath(os.curdir)

  # Clean up malformed version numbers.
  version = sanitizeVersionInput(version)

  packageName = buildPackage(sourceDir, version)
  testPackage(packageName)
  tagRepository(sourceDir, version)
# end if
