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

brickNumericdir		= $(includedir)/brick/numeric
brickNumeric_HEADERS	= \
	array1D.hh array1D_impl.hh \
	array2D.hh array2D_impl.hh \
	array3D.hh array3D_impl.hh \
	derivativeRidders.hh \
	functional.hh \
	index2D.hh \
	index3D.hh \
	maxRecorder.hh \
	minRecorder.hh \
	numericTraits.hh \
	rotations.hh rotations_impl.hh \
	slice.hh \
	solveQuadratic.hh solveQuadratic_impl.hh \
	subArray1D.hh subArray1D_impl.hh \
	subArray2D.hh subArray2D_impl.hh \
	quaternion.hh quaternion_impl.hh \
	transform2D.hh transform2D_impl.hh \
	transform3D.hh transform3D_impl.hh \
	transform3DTo2D.hh transform3DTo2D_impl.hh \
	vector2D.hh vector2D_impl.hh \
	vector3D.hh vector3D_impl.hh \
	utilities.hh utilities_impl.hh


# Define the library to build.

lib_LTLIBRARIES		= libbrickNumeric.la
libbrickNumeric_la_SOURCES = \
	index2D.cc \
	index3D.cc


# List the files that need to be cleaned.

MAINTAINERCLEANFILES    = Makefile.in
