#!/usr/bin/python

import os
import string
import sys

HOME = os.environ['HOME']
DLROPTIMIZATION = os.path.realpath(os.curdir)
GREP = '/bin/grep'
MAKE = '/usr/bin/make'
MV = '/bin/mv'
RM = '/bin/rm'
TAR = '/bin/tar'
ZIP = '/usr/bin/zip'

if __name__ == '__main__':
  if len(sys.argv) != 2:
    print 'Usage: %s version' % sys.argv[0]
    print 'Example: %s 1.11' % sys.argv[0]
    sys.exit(65)
  # end if

  VERSION = sys.argv[1]

  revision = 'version%s' % string.replace(VERSION, '.', '_')

  os.chdir(DLROPTIMIZATION)

  errorNumber = 1

  if os.system(GREP + ' "%s" VERSION.TXT' % VERSION) != 0:
    sys.exit(errorNumber)
  # end if

  if os.system('%s AC_INIT configure.ac | %s "%s"'
               % (GREP, GREP, VERSION)) != 0:
    sys.exit(errorNumber)
  # end if

  os.system(MAKE + ' maintainer-clean')

  errorNumber += 1
  if os.system('./bootstrap') != 0:
    sys.exit(errorNumber)
  # end if

  errorNumber += 1
  if os.system('./configure') != 0:
    sys.exit(errorNumber)
  # end if

  errorNumber += 1
  if os.system(MAKE + ' distcheck') != 0:
    sys.exit(errorNumber)
  # end if
  
  os.chdir('%s/tmp' % HOME)
  
  errorNumber += 1
  if os.system(TAR + ' -zxvf %s/dlroptimization-"%s".tar.gz' % (DLROPTIMIZATION, VERSION)) != 0:
    sys.exit(errorNumber)
  # end if

  os.chdir('dlroptimization-%s' % VERSION)
  
  errorNumber += 1
  if os.system('./configure') != 0:
    sys.exit(errorNumber)
  # end if
  
  errorNumber += 1
  if os.system(MAKE + ' apidoc') != 0:
    sys.exit(errorNumber)
  # end if
  
  os.chdir('doc')
  
  errorNumber += 1
  if os.system(TAR + ' -zcvf %s/dlroptimization-"%s"_htmlDoc.tgz html'
               % (DLROPTIMIZATION, VERSION)) != 0:
    sys.exit(errorNumber)
  # end if

  if 0:
    os.chdir('latex')
    
    errorNumber += 1
    if os.system(MAKE) != 0:
      sys.exit(errorNumber)
    # end if
  
    errorNumber += 1
    if os.system(MV + ' refman.pdf %s/dlroptimization-"%s"_manual.pdf'
                 % (DLROPTIMIZATION, VERSION)) != 0:
      sys.exit(errorNumber)
    # end if
  #end if
  
  os.chdir('%s/tmp' % HOME)
  
  errorNumber += 1
  if os.system(RM + ' -rf dlroptimization-"%s"' % VERSION) != 0:
    sys.exit(errorNumber)
  # end if

  os.chdir(DLROPTIMIZATION)

  errorNumber += 1
##   if os.system('cvs tag -F "%s"' % revision) != 0:
##     sys.exit(errorNumber)
##   # end if

  sys.exit(0)
#end if
