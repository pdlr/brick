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
	cameraIntrinsics.hh \
	cameraIntrinsicsDistorted.hh \
	cameraIntrinsicsPinhole.hh \
	cameraIntrinsicsPlumbBob.hh cameraIntrinsicsPlumbBob_impl.hh \
	colorspaceConverter.hh colorspaceConverter_impl.hh \
	image.hh \
	imageFilter.hh imageFilter_impl.hh \
	imageIO.hh imageIO_impl.hh \
	imageFormat.hh \
	imageFormatTraits.hh imageFormatTraits_impl.hh \
	imagePyramid.hh imagePyramid_impl.hh \
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
	utilities.hh utilities_impl.hh

stillToComeHeaders		= \
	calibrationTools.h \
	calibrationToolsRobust.h \
	canny.h \
	connectedComponents.h \
	dilate.h \
	disjointSet.h \
	eightPointAlgorithm.h \
	erode.h \
	extendedKalmanFilter.h \
	featureAssociation.h \
	fivePointAlgorithm.h \
	fourierTransformer.h \
	getEuclideanDistance.h \
	histogramEqualize.h \
	imageWarper.h \
	kdTree.h \
	nChooseKSampleSelector.h \
	naiveSnake.h \
	nonMaximumSuppress.h \
	opticalFlow.h \
	quadMap.h \
	registerPoints3D.h \
	segmenterFelzenszwalb.h \
	sobel.h \
	stereoRectify.h \
	threePointAlgorithm.h


# Define the library to build.
lib_LTLIBRARIES		= libbrickComputerVision.la
libbrickComputerVision_la_SOURCES = \
	cameraIntrinsicsPinhole.cc cameraIntrinsicsPinhole.hh \
	cameraIntrinsicsPlumbBob.cc cameraIntrinsicsPlumbBob.hh \
	cameraIntrinsicsPlumbBob_impl.hh \
	imageIO.cc imageIO.hh \
	keypointMatcherFast.cc keypointMatcherFast.hh \
	keypointSelectorFast.cc keypointSelectorFast.hh \
	ransac.cc ransac.hh

otherSources =\
	eightPointAlgorithm.cpp eightPointAlgorithm.h \
	extendedKalmanFilter.cpp extendedKalmanFilter.h \
	fivePointAlgorithm.cpp fivePointAlgorithm.h \
	fourierTransformer.cpp fourierTransformer.h \
	histogramEqualize.cpp histogramEqualize.h \
	naiveSnake.cpp naiveSnake.h \
	segmenterFelzenszwalb.cpp segmenterFelzenszwalb.h \
	stereoRectify.cpp stereoRectify.h

# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in
