#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# We need to recurse into the following directories.

SUBDIRS			= test

# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Define the list of header files to install and the location to
# install them.

brickGeometrydir		= $(includedir)/brick/geometry
brickGeometry_HEADERS	= \
	circle2D.hh circle2D_impl.hh \
	ellipse2D.hh ellipse2D_impl.hh \
	lineSegment2D.hh lineSegment2D_impl.hh \
	plane3D.hh plane3D_impl.hh \
	ray2D.hh ray2D_impl.hh \
	ray3D.hh ray3D_impl.hh \
	triangle3D.hh triangle3D_impl.hh \
	utilities2D.hh utilities3D.hh


# Define the library to build.

lib_LTLIBRARIES		= libbrickGeometry.la
libbrickGeometry_la_SOURCES = 


# List the files that need to be cleaned.

MAINTAINERCLEANFILES    = Makefile.in

