#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Specify which tests should be run.

TESTS			= \
	amanatidesWoo2DTest \
	amanatidesWoo3DTest \
	array1DTest \
	array2DTest \
	array3DTest \
	arrayNDTest \
	boxIntegrator2DTest \
	bSplineTest \
	convolve1DTest \
	convolve2DTest \
	convolveNDTest \
	derivativeRiddersTest \
	filterTest \
	ieeeFloat32Test \
	index3DTest \
	geometry2DTest \
	maxRecorderTest \
	minRecorderTest \
	normalizedCorrelatorTest \
	polynomialTest \
	rotationsTest \
	sampledFunctionsTest \
	solveCubicTest \
	solveQuadraticTest \
	solveQuarticTest \
	stencil2DTest \
	subpixelInterpolateTest \
	transform2DTest \
	transform3DTest \
	utilitiesTest


# Specify which tests should be build from source.

check_PROGRAMS          = \
	amanatidesWoo2DTest \
	amanatidesWoo3DTest \
	array1DTest \
	array2DTest \
	array3DTest \
	arrayNDTest \
	boxIntegrator2DTest \
	bSplineTest \
	convolve1DTest \
	convolve2DTest \
	convolveNDTest \
	derivativeRiddersTest \
	filterTest \
	ieeeFloat32Test \
	index3DTest \
	geometry2DTest \
	maxRecorderTest \
	minRecorderTest \
	normalizedCorrelatorTest \
	polynomialTest \
	rotationsTest \
	sampledFunctionsTest \
	solveCubicTest \
	solveQuadraticTest \
	solveQuarticTest \
	stencil2DTest \
	subpixelInterpolateTest	\
	transform2DTest \
	transform3DTest \
	utilitiesTest


# For each test which should be built from source, specify which 
# source files.

amanatidesWoo2DTest_SOURCES	= amanatidesWoo2DTest.cc
amanatidesWoo3DTest_SOURCES	= amanatidesWoo3DTest.cc
array1DTest_SOURCES		= array1DTest.cc arrayTestCommon.hh
array2DTest_SOURCES		= array2DTest.cc arrayTestCommon.hh
array3DTest_SOURCES		= array3DTest.cc arrayTestCommon.hh
arrayNDTest_SOURCES		= arrayNDTest.cc arrayTestCommon.hh
boxIntegrator2DTest_SOURCES	= boxIntegrator2DTest.cc
bSplineTest_SOURCES		= bSplineTest.cc
convolve1DTest_SOURCES		= convolve1DTest.cc
convolve2DTest_SOURCES		= convolve2DTest.cc
convolveNDTest_SOURCES		= convolveNDTest.cc
derivativeRiddersTest_SOURCES	= derivativeRiddersTest.cc
filterTest_SOURCES		= filterTest.cc
ieeeFloat32Test_SOURCES		= ieeeFloat32Test.cc
index3DTest_SOURCES		= index3DTest.cc
geometry2DTest_SOURCES		= geometry2DTest.cc
maxRecorderTest_SOURCES		= maxRecorderTest.cc
minRecorderTest_SOURCES		= minRecorderTest.cc
normalizedCorrelatorTest_SOURCES = normalizedCorrelatorTest.cc
polynomialTest_SOURCES		= polynomialTest.cc
rotationsTest_SOURCES		= rotationsTest.cc
sampledFunctionsTest_SOURCES	= sampledFunctionsTest.cc
solveCubicTest_SOURCES		= solveCubicTest.cc
solveQuadraticTest_SOURCES	= solveQuadraticTest.cc
solveQuarticTest_SOURCES	= solveQuarticTest.cc
stencil2DTest_SOURCES		= stencil2DTest.cc
subpixelInterpolateTest_SOURCES	= subpixelInterpolateTest.cc
transform2DTest_SOURCES		= transform2DTest.cc arrayTestCommon.hh
transform3DTest_SOURCES		= transform3DTest.cc arrayTestCommon.hh
utilitiesTest_SOURCES		= utilitiesTest.cc


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Explicitly set library locations to avoid linking with previously
# installed versions of the libraries.

AM_LDFLAGS		= -L$(top_builddir)/brick/numeric


# Specify which libraries are needed to build these tests.

LDADD			= -lbrickNumeric \
			  -lbrickPortability \
			  -lbrickCommon \
			  -lbrickTest \
			  -lbrickTestAutoMain


# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in
