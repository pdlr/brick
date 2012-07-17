#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# We need to recurse into the following directories.

SUBDIRS = test


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Define the list of header files to install and the location to
# install them.

brickComputerVisiondir		= $(includedir)/brick/computerVision
brickComputerVision_HEADERS	= \
	calibrationTools.hh calibrationTools_impl.hh \
	calibrationToolsRobust.hh calibrationToolsRobust_impl.hh \
	cameraIntrinsics.hh \
	cameraIntrinsicsDistorted.hh \
	cameraIntrinsicsPinhole.hh cameraIntrinsicsPinhole_impl.hh \
	cameraIntrinsicsPlumbBob.hh cameraIntrinsicsPlumbBob_impl.hh \
	canny.hh canny_impl.hh \
	colorspaceConverter.hh colorspaceConverter_impl.hh \
	connectedComponents.hh connectedComponents_impl.hh \
	dilate.hh dilate_impl.hh \
	disjointSet.hh disjointSet_impl.hh \
	eightPointAlgorithm.hh eightPointAlgorithm_impl.hh \
	erode.hh erode_impl.hh \
	extendedKalmanFilter.hh extendedKalmanFilter_impl.hh \
	featureAssociation.hh featureAssociation_impl.hh \
	fivePointAlgorithm.hh fivePointAlgorithm_impl.hh \
	getEuclideanDistance.hh getEuclideanDistance_impl.hh \
	histogramEqualize.hh \
	image.hh \
	imageFilter.hh imageFilter_impl.hh \
	imageIO.hh imageIO_impl.hh \
	imageFormat.hh \
	imageFormatTraits.hh imageFormatTraits_impl.hh \
	imagePyramid.hh imagePyramid_impl.hh \
	imageWarper.hh imageWarper_impl.hh \
	iterativeClosestPoint.hh iterativeClosestPoint_impl.hh \
	kdTree.hh kdTree_impl.hh \
	kernel.hh kernel_impl.hh \
	kernels.hh kernels_impl.hh \
	keypointMatcherFast.hh keypointMatcherFast_impl.hh \
	keypointSelectorFast.hh keypointSelectorFast_impl.hh \
	keypointSelectorHarris.hh keypointSelectorHarris_impl.hh \
	naiveSnake.hh naiveSnake_impl.hh \
	nonMaximumSuppress.hh nonMaximumSuppress_impl.hh \
	nChooseKSampleSelector.hh nChooseKSampleSelector_impl.hh \
	pixelBGRA.hh \
	pixelHSV.hh \
	pixelRGB.hh \
	pixelRGBA.hh \
	pixelYIQ.hh \
	randomSampleSelector.hh randomSampleSelector_impl.hh \
	ransac.hh ransac_impl.hh \
	ransacClassInterface.hh \
	registerPoints3D.hh registerPoints3D_impl.hh \
	segmenterFelzenszwalb.hh segmenterFelzenszwalb_impl.hh \
	sobel.hh sobel_impl.hh \
	stereoRectify.hh stereoRectify_impl.hh \
	threePointAlgorithm.hh threePointAlgorithm_impl.hh \
	utilities.hh utilities_impl.hh

stillToComeHeaders		= \
	opticalFlow.h \
	quadMap.h


# Define the library to build.
lib_LTLIBRARIES		= libbrickComputerVision.la
libbrickComputerVision_la_SOURCES = \
	cameraIntrinsicsPlumbBob_impl.hh \
	imageIO.cc imageIO.hh \
	histogramEqualize.cc histogramEqualize.hh \
	keypointMatcherFast.cc keypointMatcherFast.hh \
	keypointSelectorFast.cc keypointSelectorFast.hh \
	ransac.cc ransac.hh


# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in
