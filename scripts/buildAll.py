#!/usr/bin/python

import glob
import os
import sys

# These libraries are arranged in order of dependency.  Each library
# depends only on those libraries that preceed it in the list.
libraryList = ['brickCommon',
               'brickTest',
               'brickPortability',
               'brickRandom',
               'brickNumeric',
               'brickLinearAlgebra',
               'brickOptimization',
               'brickUtilities',
               'brickThread',
               'brickGeometry',
               'brickComputerVision',
               'brickMedia']


preferredConfig = './configure'

def checkVersionEqual(version0, version1):
  if ((version0[0] == version1[0]) & (version0[1] == version1[1]) & (version0[2] == version1[2])):
    return True
  # end if
  return False
# end if


def checkVersionGreater(version0, version1):
  if version0[0] > version1[0]:
    return True
  # end if
  if ((version0[0] == version1[0]) & (version0[1] > version1[1])):
    return True
  # end if
  if ((version0[0] == version1[0]) & (version0[1] == version1[1])
      & (version0[2] > version1[2])):
    return True
  # end if
  return False
# end if


def findLatestVersion(fileNameList):
  latestVersion = (0, 0, '')
  latestFileName = ''
  for fileName in fileNameList:
    candidateName = fileName
    while 1:
      base, ext = os.path.splitext(candidateName)
      if ((ext != '.gz') & (ext != '.tar') & (ext != '.tgz')):
        break
      # end if
      candidateName = base
    # end while
    if len(ext) <= 1:
      raise IOError('Malformed fileName: %s' % fileName)
    # end if
    try:
      if (len(ext) == 2):
        minorRev = int(ext[1:])
        revModifier = ''
      else:
        # Warning(xxx): assuming minor revision is only single digit.
        minorRev = int(ext[1:2])
        revModifier = ext[2:]
      # end if
      ext = base.split('-')[-1]
      majorRev = int(ext)
    except ValueError:
      raise IOError('Malformed fileName: %s' % fileName)
    # end try
    newVersion = (majorRev, minorRev, revModifier)
    if checkVersionGreater(newVersion, latestVersion):
      latestVersion = newVersion
      latestFileName = fileName
    # end if
  # end for
  return latestFileName, latestVersion
# end def
      

def findSource(libraryName):
  sourceDirGlob = libraryName.lower() + '-?.*'
  tarballGlob = sourceDirGlob + '.tar.gz'

  sourceDirList = glob.glob(sourceDirGlob)
  tarballList = glob.glob(tarballGlob)

  tarballName, tarballVersion = findLatestVersion(tarballList)
  prevSourceDirName, prevSourceVersion = findLatestVersion(sourceDirList)

  isTarballExisting = (tarballName != '')
  isSourceDirExisting = (prevSourceDirName != '')

  if isTarballExisting:
    sourceDirName = os.path.splitext(os.path.splitext(tarballName)[0])[0]
  else:
    sourceDirName = prevSourceDirName
  # end if

  if (isTarballExisting & isSourceDirExisting):
    if checkVersionGreater(prevSourceVersion, tarballVersion):
      raise IOException(
        'Tarball %s is out of date relative to source dir %s'
        % (tarballName, prevSourceDirName))
    # end if
    if not checkVersionEqual(prevSourceVersion, tarballVersion):
      # Source dir is out of date.
      isSourceDirExisting = False
    # end if
  # end if

  return (tarballName, sourceDirName, isTarballExisting, isSourceDirExisting)
# end def


def unpackTarball(tarballName, sourceDirName):
  if tarballName != '':
    if os.path.exists(sourceDirName):
      print 'Source directory:\n  %s\nalready exists.\n' % sourceDirName
      response = raw_input('Overwrite? [y/N]: ')
      if len(response) == 0:
        response = 'N'
      # end if
      if ((response[0] == 'N') | (response[0] == 'n')):
        return False
      # end if
      returnCode = os.system('rm -rf %s' % sourceDirName)
      if returnCode != 0:
        raise IOError('Couldn\'t remove old source tree: %s' % sourceDirName)
      # end if
    # end if
    returnCode = os.system('tar -zxf %s' % tarballName)
    if returnCode != 0:
      raise IOError('Couldn\'t unpack tarball: %s' % tarballName)
    # end if
  # end if
  return True
# end def


if __name__ == '__main__':
  if len(sys.argv) < 2:
    print "Usage: %s [--clean] installPrefix" % sys.argv[0]
    print ""
    print "  where installPrefix is the directory under which you'd like"
    print "  all of the libraries, headers, etc. to be installed."
    sys.exit(0)
  # end if

  installPrefix = "."
  action = 'build'
  argIndex = 1
  while argIndex < len(sys.argv):
    if sys.argv[argIndex][0] == '-':
      if sys.argv[argIndex] == '--clean':
        action = 'clean'
      else:
        raise ValueError("Unrecognized option: %s" % sys.argv[argIndex])
      # end if
    else:
      installPrefix = sys.argv[argIndex]
    # end if
    argIndex = argIndex + 1
  # end if

  installPrefix = os.path.abspath(installPrefix)

  for libraryName in libraryList:
    print '============= %s ============= ' % libraryName
    if os.path.isdir(libraryName):
      print '  Found VCS checkout for %s...' % libraryName
      response = raw_input('  Rebuild? [Y/n]: ')
      if response == '':
        response = 'Y'
      # end if
      if ((response[0] == 'y') | (response[0] == 'Y')):
        os.system('cd %s && make maintainer-clean' % (libraryName))
        returnCode = os.system(
          'cd %s && ./bootstrap && sh %s && make && make dist && make install' 
          % (libraryName, preferredConfig))
        if returnCode != 0:
          raise ValueError(
            'Couldn\'t build %s from VCS checkout.' % libraryName)
        # end if
      # end if
    else:
      # Assume we're going to compile!
      isOKToContinue = True

      # Find out what materials we have to work with.
      tarballName, sourceDirName, isTarballExisting, isSourceDirExisting = (
        findSource(libraryName))

      # If there's a tarball, unpack it, unless the user tells you not
      # to overwrite an existing source tree, in which case do
      # nothing.
      if isTarballExisting:
        isOKToContinue = unpackTarball(tarballName, sourceDirName)
        isSourceDirExisting = True
      # end if

      # Do the actual build.
      if (isSourceDirExisting & isOKToContinue):
        cflags = "-I%s/include" % installPrefix
        cxxflags = "-I%s/include" % installPrefix
        if not os.environ.has_key("CFLAGS"):
          cflags = cflags + " -O2"
        # end if
        if not os.environ.has_key("CXXFLAGS"):
          cxxflags = cxxflags + " -O2"
        # end if
        returnCode = os.system(
          """export LDFLAGS="$LDFLAGS -L%s/lib" && \
             export CFLAGS="$CFLAGS %s" && \
             export CXXFLAGS="$CXXFLAGS %s" && \
             cd %s && \
             ./configure --prefix "%s" && \
             make && \
             make check && \
             make install""" 
          % (installPrefix, cflags, cxxflags, sourceDirName, installPrefix))
        if returnCode != 0:
          # Clean up failed build.
          print 'Failed build...'
          cleanupCommand = 'rm -rf %s' % sourceDirName
          response = raw_input('  Clean up (%s)? [y/N]: ' % cleanupCommand)
          if response == '':
            response = 'N'
          # end if
          if ((response[0] == 'y') | (response[0] == 'Y')):
            os.system(cleanupCommand)
          # end if
          raise ValueError(
            'Couldn\'t build %s from tarball or existing source.'
            % libraryName)
        # end if
      # end if
    # end if
  # end for
# end if
