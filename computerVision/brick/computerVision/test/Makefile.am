#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#
# $Revision: 1335 $
# $Date: 2010-07-03 14:20:37 -0400 (Sat, 03 Jul 2010) $
#


# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Additional files which should be included in the distribution.

EXTRA_DIST              = \
	connectedComponentsTestImage0.pgm \
	dilateErodeTestImage0.pgm \
	dilatedTestImage0.pgm \
	edgeTestImage0.pgm \
	erodedTestImage0.pgm \
	testImagePGM0.pgm \
	testImagePGM1.pgm \
	testImagePPM0.ppm \
	testImagePPM1.ppm


# Specify which tests should be run.

TESTS = \
	colorspaceConverterTest \
	imageFilterTest \
	imagePyramidTest \
	keypointSelectorLepetitTest \
	ransacTest \
	utilitiesTest


# Specify which tests should be build from source.

check_PROGRAMS = \
	colorspaceConverterTest \
	imageFilterTest \
	imagePyramidTest \
	keypointSelectorLepetitTest \
	ransacTest \
	utilitiesTest


# For each test which should be built from source, specify which 
# source files.

colorspaceConverterTest_SOURCES	= colorspaceConverterTest.cc
imageFilterTest_SOURCES		= imageFilterTest.cc
imagePyramidTest_SOURCES	= imagePyramidTest.cc
keypointSelectorLepetitTest_SOURCES	= keypointSelectorLepetitTest.cc
ransacTest_SOURCES		= ransacTest.cc
utilitiesTest_SOURCES		= utilitiesTest.cc


# Set the include path and specify where to look for data files.

AM_CPPFLAGS		= -I$(top_srcdir) -DBRICK_TEST_DATA_DIR="\"$(srcdir)\""


# Explicitly set library locations to avoid linking with previously
# installed versions of the libraries.

AM_LDFLAGS		= -L$(top_builddir)/brick/computerVision


# Specify which libraries are needed to build these tests.

LDADD			= -lbrickComputerVision \
			  -lbrickUtilities \
			  -lbrickGeometry \
			  -lbrickRandom \
			  -lbrickLinearAlgebra \
			  -lbrickNumeric \
			  -lbrickPortability \
			  -lbrickCommon \
			  -lbrickTest \
			  -lbrickTestAutoMain \
			  -llapack -lblas \
			  -lpng \
			  -lfftw3 \
			  -lm



# List the files that need to be cleaned.

CLEANFILES		= foo.pgm foo.png
MAINTAINERCLEANFILES	= Makefile.in
