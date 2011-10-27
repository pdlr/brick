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
	randomSampleSelector.hh randomSampleSelector_impl.hh\
	ransac.hh ransac_impl.hh

stillToComeHeaders		= \
	calibrationTools.h \
	calibrationToolsRobust.h \
	cameraIntrinsics.h \
	cameraIntrinsicsDistorted.h \
	cameraIntrinsicsPinhole.h \
	cameraIntrinsicsPlumbBob.h \
	canny.h \
	colorspaceConverter.h \
	connectedComponents.h \
	dilate.h \
	disjointSet.h \
	eightPointAlgorithm.h \
	erode.h \
	extendedKalmanFilter.h \
	featureAssociation.h \
	filter.h \
	fivePointAlgorithm.h \
	fourierTransformer.h \
	getEuclideanDistance.h \
	histogramEqualize.h \
	image.h \
	imageFormat.h \
	imageFormatTraits.h \
	imageIO.h \
	imageWarper.h \
	kdTree.h \
	kernel.h \
	kernels.h \
	nChooseKSampleSelector.h \
	naiveSnake.h \
	nonMaximumSuppress.h \
	opticalFlow.h \
	pixelBGRA.h \
	pixelHSV.h \
	pixelRGB.h \
	pixelRGBA.h \
	pixelYIQ.h \
	quadMap.h \
	registerPoints3D.h \
	segmenterFelzenszwalb.h \
	sobel.h \
	stereoRectify.h \
	threePointAlgorithm.h \
	utilities.h


# Define the library to build.
lib_LTLIBRARIES		= libbrickComputerVision.la
libbrickComputerVision_la_SOURCES = 

otherSources =\
	cameraIntrinsicsPinhole.cpp cameraIntrinsicsPinhole.h \
	cameraIntrinsicsPlumbBob.cpp cameraIntrinsicsPlumbBob.h \
	eightPointAlgorithm.cpp eightPointAlgorithm.h \
	extendedKalmanFilter.cpp extendedKalmanFilter.h \
	fivePointAlgorithm.cpp fivePointAlgorithm.h \
	fourierTransformer.cpp fourierTransformer.h \
	histogramEqualize.cpp histogramEqualize.h \
	imageIO.cpp imageIO.h \
	naiveSnake.cpp naiveSnake.h \
	segmenterFelzenszwalb.cpp segmenterFelzenszwalb.h \
	stereoRectify.cpp stereoRectify.h

# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in
