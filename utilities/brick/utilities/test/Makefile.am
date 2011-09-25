#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Specify which tests should be run.

TESTS			= \
	dateTest \
	imageIOTest \
	lockFileTest \
	optionParserTest \
	pathTest \
	stringManipulationTest \
	timeUtilitiesTest


# Specify which tests should be build from source.

check_PROGRAMS          = \
	dateTest \
	imageIOTest \
	lockFileTest \
	optionParserTest \
	pathTest \
	stringManipulationTest \
	timeUtilitiesTest


# For each test which should be built from source, specify which 
# source files.

dateTest_SOURCES		= dateTest.cc
imageIOTest_SOURCES		= imageIOTest.cc
lockFileTest_SOURCES		= lockFileTest.cc
optionParserTest_SOURCES	= optionParserTest.cc
pathTest_SOURCES		= pathTest.cc
stringManipulationTest_SOURCES	= stringManipulationTest.cc
timeUtilitiesTest_SOURCES	= timeUtilitiesTest.cc


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir) -DBRICK_TEST_USE_AUTOMATIC_MAIN


# Explicitly set library locations to avoid linking with previously
# installed versions of the libraries.

AM_LDFLAGS		= -L$(top_builddir)/brick/utilities


# Specify which libraries are needed to build these tests.

LDADD			= \
	-lbrickUtilities \
	-lbrickCommon \
	-lbrickPortability \
	-lbrickTest \
	-lbrickTestAutoMain


# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in bar.pgm bar.ppm foo.pgm foo.ppm
