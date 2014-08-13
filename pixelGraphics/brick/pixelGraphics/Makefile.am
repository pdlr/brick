#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# Set the include path.

AM_CPPFLAGS		= -I$(top_srcdir)


# Define the list of header files to install and the location to
# install them.

brickPixelGraphicsdir	= $(includedir)/brick/pixelGraphics
brickPixelGraphics_HEADERS = \
	draw2D.hh draw2D_impl.hh


# Define the library to build.

#   No library to build, everything is contained in the header files.


# List the files that need to be cleaned.

MAINTAINERCLEANFILES    = Makefile.in
