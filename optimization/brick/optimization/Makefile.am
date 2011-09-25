#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Define the list of header files to install and the location to
# install them.

brickOptimizationdir	= $(includedir)/brickOptimization
brickOptimization_HEADERS = \
	gradientFunction.h \
	gradientFunctionLM.h \
	optimizer.h \
	optimizerBFGS.h \
	optimizerCommon.h \
	optimizerLM.h \
	optimizerLineSearch.h \
	optimizerNelderMead.h


# Define the library to build.

#   No library to build, everything is contained in the header files.


# List the files that need to be cleaned.

MAINTAINERCLEANFILES    = Makefile.in

