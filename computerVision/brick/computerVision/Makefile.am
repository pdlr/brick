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
	colorspaceConverter.hh colorspaceConverter_impl.hh \
	connectedComponents.hh connectedComponents_impl.hh \
	dilate.hh dilate_impl.hh \
	disjointSet.hh disjointSet_impl.hh \
	eightPointAlgorithm.hh eightPointAlgorithm_impl.hh \
	erode.hh erode_impl.hh \
	extendedKalmanFilter.hh extendedKalmanFilter_impl.hh \
	fivePointAlgorithm.hh fivePointAlgorithm_impl.hh \
	image.hh \
	imageFilter.hh imageFilter_impl.hh \
	imageIO.hh imageIO_impl.hh \
	imageFormat.hh \
	imageFormatTraits.hh imageFormatTraits_impl.hh \
	imagePyramid.hh imagePyramid_impl.hh \
	imageWarper.hh imageWarper_impl.hh \
	kdTree.hh kdTree_impl.hh \
	kernel.hh kernel_impl.hh \
	kernels.hh kernels_impl.hh \
	keypointMatcherFast.hh keypointMatcherFast_impl.hh \
	keypointSelectorFast.hh keypointSelectorFast_impl.hh \
	pixelBGRA.hh \
	pixelHSV.hh \
	pixelRGB.hh \
	pixelRGBA.hh \
	pixelYIQ.hh \
	randomSampleSelector.hh randomSampleSelector_impl.hh \
	ransac.hh ransac_impl.hh \
	registerPoints3D.hh registerPoints3D_impl.hh \
	segmenterFelzenszwalb.hh segmenterFelzenszwalb_impl.hh \
	stereoRectify.hh stereoRectify_impl.hh \
	threePointAlgorithm.hh threePointAlgorithm_impl.hh \
	utilities.hh utilities_impl.hh

stillToComeHeaders		= \
	canny.h \
	featureAssociation.h \
	fourierTransformer.h \
	getEuclideanDistance.h \
	histogramEqualize.h \
	nChooseKSampleSelector.h \
	naiveSnake.h \
	nonMaximumSuppress.h \
	opticalFlow.h \
	quadMap.h \
	sobel.h


# Define the library to build.
lib_LTLIBRARIES		= libbrickComputerVision.la
libbrickComputerVision_la_SOURCES = \
	cameraIntrinsicsPlumbBob_impl.hh \
	imageIO.cc imageIO.hh \
	keypointMatcherFast.cc keypointMatcherFast.hh \
	keypointSelectorFast.cc keypointSelectorFast.hh \
	ransac.cc ransac.hh

otherSources =\
	extendedKalmanFilter.cpp extendedKalmanFilter.h \
	fourierTransformer.cpp fourierTransformer.h \
	histogramEqualize.cpp histogramEqualize.h \
	naiveSnake.cpp naiveSnake.h \
	stereoRectify.cpp stereoRectify.h

# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in
