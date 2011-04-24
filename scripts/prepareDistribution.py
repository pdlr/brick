#!/usr/bin/python

import os
import string
import sys

GREP = '/bin/grep'
MAKE = '/usr/bin/make'
MV = '/bin/mv'
RM = '/bin/rm'
TAR = '/bin/tar'
ZIP = '/usr/bin/zip'
TEMPDIR = '/var/tmp/brickBuild'

def prepareDistribution(packageName, sourceDir, version):

  # Tarfiles are named "brickfoo.tar.gz", not "brickFoo.tar.gz", so we
  # need a lowercase version of the package name.
  packageNameLC = packageName.lower()

  # We need a civilized version of the revision number to use as a tag
  # in the VCS.
  revisionTag = ('%s_version%s' % 
                 (packageName, string.replace(version, '.', '_')))

  # On with the show.
  os.chdir(sourceDir)

  errorNumber = 1

  if os.system(GREP + ' "%s" VERSION.TXT' % version) != 0:
    return errorNumber
  # end if

  errorNumber += 1
  if os.system('%s AC_INIT configure.ac | %s "%s"'
               % (GREP, GREP, version)) != 0:
    return errorNumber
  # end if

  os.system(MAKE + ' maintainer-clean')

  errorNumber += 1
  if os.system('./bootstrap') != 0:
    return errorNumber
  # end if

  errorNumber += 1
  if os.system('./configure') != 0:
    return errorNumber
  # end if

  errorNumber += 1
  if os.system(MAKE + ' distcheck') != 0:
    return errorNumber
  # end if
  
  if not os.path.isdir(TEMPDIR):
    if os.path.exists(TEMPDIR):
      os.system('rm -f "%s"' % TEMPDIR)
    # end if
    os.system('mkdir -p "%s"' % TEMPDIR)
  # end if
  os.chdir(TEMPDIR)
  
  errorNumber += 1
  if os.system(TAR + ' -zxvf %s/%s-"%s".tar.gz' 
               % (sourceDir, packageNameLC, version)) != 0:
    return errorNumber
  # end if

  os.chdir('%s-%s' % (packageNameLC, version))
  
  errorNumber += 1
  if os.system('./configure') != 0:
    return errorNumber
  # end if
  
  errorNumber += 1
  if os.system(MAKE + ' apidoc') != 0:
    return errorNumber
  # end if
  
  os.chdir('doc')
  
  errorNumber += 1
  if os.system(TAR + ' -zcvf %s/%s-"%s"_htmlDoc.tgz html'
               % (sourceDir, packageNameLC, version)) != 0:
    return errorNumber
  # end if

  if 0:
    os.chdir('latex')
    
    errorNumber += 1
    if os.system(MAKE) != 0:
      return errorNumber
    # end if
  
    errorNumber += 1
    if os.system(MV + ' refman.pdf %s/%s-"%s"_manual.pdf'
                 % (sourceDir, packageNameLC, version)) != 0:
      return errorNumber
    # end if
  #end if
  
  os.chdir(TEMPDIR)
  
  errorNumber += 1
  if os.system(RM + ' -rf %s-"%s"' % (packageNameLC, version)) != 0:
    return errorNumber
  # end if

  os.chdir(sourceDir)

  errorNumber += 1
  if os.system('hg tag -f "%s"' % revisionTag) != 0:
    return errorNumber
  # end if

  return 0
# end def


if __name__ == '__main__':
  if len(sys.argv) != 3:
    print 'Usage: %s packageName version' % sys.argv[0]
    print 'Example: %s brickCommon 1.11' % sys.argv[0]
    sys.exit(65)
  # end if

  packageName = sys.argv[1]
  version = sys.argv[2]
  sourceDir = os.path.realpath(os.curdir)
  
  returnCode = prepareDistribution(packageName, sourceDir, version)
  sys.exit(returnCode)
# end if
