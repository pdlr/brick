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

EXTRA_DIST              = connectedComponentsTestImage0.pgm \
			  dilateErodeTestImage0.pgm \
			  dilatedTestImage0.pgm \
			  edgeTestImage0.pgm \
			  erodedTestImage0.pgm \
			  testImagePGM0.pgm \
			  testImagePGM1.pgm \
			  testImagePGM2.pgm \
			  testImagePGM3.pgm \
			  testImagePPM0.ppm \
			  testImagePPM1.ppm


# Specify which tests should be run.

TESTS			= calibrationToolsTest \
			  cameraIntrinsicsPinholeTest \
			  cameraIntrinsicsPlumbBobTest \
			  cannyTest \
			  colorspaceConverterTest \
			  connectedComponentsTest \
			  dilateTest \
			  eightPointAlgorithmTest \
			  erodeTest \
			  extendedKalmanFilterTest \
			  featureAssociationTest \
			  filterTest \
			  fivePointAlgorithmTest \
			  fourierTransformerTest \
			  getEuclideanDistanceTest \
			  imageIOTest \
			  imageWarperTest \
			  kdTreeTest \
			  nChooseKSampleSelectorTest \
			  naiveSnakeTest \
			  nonMaximumSuppressTest \
			  opticalFlowTest \
			  quadMapTest \
			  ransacTest \
			  registerPoints3DTest \
			  segmenterFelzenszwalbTest \
			  sobelTest \
			  stereoRectifyTest \
			  threePointAlgorithmTest \
			  utilitiesTest


# Specify which tests should be build from source.

check_PROGRAMS          = calibrationToolsTest \
			  cameraIntrinsicsPinholeTest \
			  cameraIntrinsicsPlumbBobTest \
			  cannyTest \
			  colorspaceConverterTest \
			  connectedComponentsTest \
			  dilateTest \
			  eightPointAlgorithmTest \
			  erodeTest \
			  extendedKalmanFilterTest \
			  featureAssociationTest \
			  filterTest \
			  fivePointAlgorithmTest \
			  fourierTransformerTest \
			  getEuclideanDistanceTest \
			  imageIOTest \
			  imageWarperTest \
			  kdTreeTest \
			  nChooseKSampleSelectorTest \
			  naiveSnakeTest \
			  nonMaximumSuppressTest \
			  opticalFlowTest \
			  quadMapTest \
			  ransacTest \
			  registerPoints3DTest \
			  segmenterFelzenszwalbTest \
			  sobelTest \
			  stereoRectifyTest \
			  threePointAlgorithmTest \
			  utilitiesTest


# For each test which should be built from source, specify which 
# source files.

calibrationToolsTest_SOURCES	= calibrationToolsTest.cpp
cameraIntrinsicsPinholeTest_SOURCES	= cameraIntrinsicsPinholeTest.cpp
cameraIntrinsicsPlumbBobTest_SOURCES	= cameraIntrinsicsPlumbBobTest.cpp
cannyTest_SOURCES		= cannyTest.cpp testImages.h
colorspaceConverterTest_SOURCES	= colorspaceConverterTest.cpp
connectedComponentsTest_SOURCES	= connectedComponentsTest.cpp
dilateTest_SOURCES		= dilateTest.cpp
eightPointAlgorithmTest_SOURCES	= eightPointAlgorithmTest.cpp
erodeTest_SOURCES		= erodeTest.cpp
extendedKalmanFilterTest_SOURCES	= extendedKalmanFilterTest.cpp
featureAssociationTest_SOURCES	= featureAssociationTest.cpp
filterTest_SOURCES		= filterTest.cpp testImages.h
fivePointAlgorithmTest_SOURCES	= fivePointAlgorithmTest.cpp
fourierTransformerTest_SOURCES	= fourierTransformerTest.cpp
getEuclideanDistanceTest_SOURCES	= getEuclideanDistanceTest.cpp
imageIOTest_SOURCES		= imageIOTest.cpp testImages.h
imageWarperTest_SOURCES		= imageWarperTest.cpp
kdTreeTest_SOURCES		= kdTreeTest.cpp
// mutualInformationTest_SOURCES	= mutualInformationTest.cpp
nChooseKSampleSelectorTest_SOURCES	= nChooseKSampleSelectorTest.cpp
naiveSnakeTest_SOURCES		= naiveSnakeTest.cpp
nonMaximumSuppressTest_SOURCES	= nonMaximumSuppressTest.cpp 
opticalFlowTest_SOURCES		= opticalFlowTest.cpp
quadMapTest_SOURCES		= quadMapTest.cpp
ransacTest_SOURCES		= ransacTest.cpp
registerPoints3DTest_SOURCES	= registerPoints3DTest.cpp
segmenterFelzenszwalbTest_SOURCES = segmenterFelzenszwalbTest.cpp
sobelTest_SOURCES		= sobelTest.cpp
stereoRectifyTest_SOURCES	= stereoRectifyTest.cpp
threePointAlgorithmTest_SOURCES	= threePointAlgorithmTest.cpp
utilitiesTest_SOURCES		= utilitiesTest.cpp testImages.h


# Set the include path and specify where to look for data files.

AM_CPPFLAGS		= -I$(top_srcdir) -DDLR_TEST_DATA_DIR="\"$(srcdir)\""


# Explicitly set library locations to avoid linking with previously
# installed versions of the libraries.

AM_LDFLAGS		= -L$(top_builddir)/dlrComputerVision


# Specify which libraries are needed to build these tests.

LDADD			= -ldlrComputerVision \
			  -ldlrUtilities \
			  -ldlrGeometry \
			  -ldlrRandom \
			  -ldlrLinearAlgebra \
			  -ldlrNumeric \
			  -ldlrPortability \
			  -ldlrCommon \
			  -ldlrTest \
			  -ldlrTestAutoMain \
			  -llapack -lblas \
			  -lpng \
			  -lfftw3 \
			  -lm



# List the files that need to be cleaned.

CLEANFILES		= foo.pgm foo.png
MAINTAINERCLEANFILES	= Makefile.in
