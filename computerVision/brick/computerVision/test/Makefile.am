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
	bullseyeTestImage0.pgm \
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
	calibrationToolsTest \
	calibrationToolsRobustTest \
	cameraIntrinsicsPinholeTest \
	cameraIntrinsicsPlumbBobTest \
	cameraIntrinsicsRationalTest \
	cannyTest \
	colorspaceConverterTest \
	connectedComponentsTest \
	dilateTest \
	eightPointAlgorithmTest \
	erodeTest \
	extendedKalmanFilterTest \
	featureAssociationTest \
	fivePointAlgorithmTest \
	getEuclideanDistanceTest \
	imageFilterTest \
	imageIOTest \
	imagePyramidTest \
	imagePyramidBinomialTest \
	imageWarperTest \
	kdTreeTest \
	keypointMatcherFastTest \
	keypointSelectorBullseyeTest \
	keypointSelectorFastTest \
	keypointSelectorHarrisTest \
	naiveSnakeTest \
	nChooseKSampleSelectorTest \
	nonMaximumSuppressTest \
	registerPoints3DTest \
	segmenterFelzenszwalbTest \
	sobelTest \
	stereoRectifyTest \
	threePointAlgorithmTest \
	utilitiesTest


# Specify which tests should be build from source.

check_PROGRAMS = \
	calibrationToolsTest \
	calibrationToolsRobustTest \
	cameraIntrinsicsPinholeTest \
	cameraIntrinsicsPlumbBobTest \
	cameraIntrinsicsRationalTest \
	cannyTest \
	colorspaceConverterTest \
	connectedComponentsTest \
	dilateTest \
	eightPointAlgorithmTest \
	erodeTest \
	extendedKalmanFilterTest \
	featureAssociationTest \
	fivePointAlgorithmTest \
	getEuclideanDistanceTest \
	imageFilterTest \
	imageIOTest \
	imagePyramidTest \
	imagePyramidBinomialTest \
	imageWarperTest \
	kdTreeTest \
	keypointMatcherFastTest \
	keypointSelectorBullseyeTest \
	keypointSelectorFastTest \
	keypointSelectorHarrisTest \
	naiveSnakeTest \
	nChooseKSampleSelectorTest \
	nonMaximumSuppressTest \
	registerPoints3DTest \
	segmenterFelzenszwalbTest \
	sobelTest \
	stereoRectifyTest \
	threePointAlgorithmTest \
	utilitiesTest


# For each test which should be built from source, specify which 
# source files.

calibrationToolsTest_SOURCES	= calibrationToolsTest.cc
calibrationToolsRobustTest_SOURCES	= calibrationToolsRobustTest.cc
cameraIntrinsicsPinholeTest_SOURCES	= cameraIntrinsicsPinholeTest.cc
cameraIntrinsicsPlumbBobTest_SOURCES	= cameraIntrinsicsPlumbBobTest.cc
cameraIntrinsicsRationalTest_SOURCES	= cameraIntrinsicsRationalTest.cc
cannyTest_SOURCES		= cannyTest.cc
colorspaceConverterTest_SOURCES	= colorspaceConverterTest.cc
connectedComponentsTest_SOURCES	= connectedComponentsTest.cc testImages.hh
dilateTest_SOURCES		= dilateTest.cc testImages.hh
eightPointAlgorithmTest_SOURCES	= eightPointAlgorithmTest.cc
erodeTest_SOURCES		= erodeTest.cc testImages.hh
extendedKalmanFilterTest_SOURCES	= extendedKalmanFilterTest.cc
featureAssociationTest_SOURCES	= featureAssociationTest.cc
fivePointAlgorithmTest_SOURCES	= fivePointAlgorithmTest.cc
getEuclideanDistanceTest_SOURCES	= getEuclideanDistanceTest.cc
imageFilterTest_SOURCES		= imageFilterTest.cc testImages.hh
imageIOTest_SOURCES		= imageIOTest.cc testImages.hh
imagePyramidTest_SOURCES	= imagePyramidTest.cc testImages.hh
imagePyramidBinomialTest_SOURCES	= imagePyramidBinomialTest.cc testImages.hh
imageWarperTest_SOURCES		= imageWarperTest.cc
kdTreeTest_SOURCES		= kdTreeTest.cc
keypointMatcherFastTest_SOURCES	= keypointMatcherFastTest.cc
keypointSelectorBullseyeTest_SOURCES	= keypointSelectorBullseyeTest.cc
keypointSelectorFastTest_SOURCES	= keypointSelectorFastTest.cc
keypointSelectorHarrisTest_SOURCES	= keypointSelectorHarrisTest.cc
naiveSnakeTest_SOURCES		= naiveSnakeTest.cc
nChooseKSampleSelectorTest_SOURCES	= nChooseKSampleSelectorTest.cc
nonMaximumSuppressTest_SOURCES	= nonMaximumSuppressTest.cc
# ransacTest_SOURCES		= ransacTest.cc
registerPoints3DTest_SOURCES	= registerPoints3DTest.cc
segmenterFelzenszwalbTest_SOURCES	= segmenterFelzenszwalbTest.cc testImages.hh
sobelTest_SOURCES		= sobelTest.cc
stereoRectifyTest_SOURCES	= stereoRectifyTest.cc
threePointAlgorithmTest_SOURCES	= threePointAlgorithmTest.cc
utilitiesTest_SOURCES		= utilitiesTest.cc testImages.hh


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
			  -lm



# List the files that need to be cleaned.

CLEANFILES		= foo.pgm foo.png
MAINTAINERCLEANFILES	= Makefile.in
