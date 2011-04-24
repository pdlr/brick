#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#
# $Revision: 915 $
# $Date: 2007-05-18 01:27:47 -0400 (Fri, 18 May 2007) $
#


AUTOMAKE_OPTIONS	= dist-zip
ACLOCAL_AMFLAGS		= -I m4

# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


# We need to recurse into the following directories.

SUBDIRS			= doc brick/utilities


# Specify the list of additional files that need to be included in the
# distribution.

AUX_DIST		= $(ac_aux_dir)/config.guess \
			  $(ac_aux_dir)/config.sub \
			  $(ac_aux_dir)/depcomp \
			  $(ac_aux_dir)/install-sh \
			  $(ac_aux_dir)/ltmain.sh \
			  $(ac_aux_dir)/missing \
			  $(ac_aux_dir)/mkinstalldirs

EXTRA_DIST              = LICENSE.TXT \
			  README.TXT \
			  RELEASE_NOTES.TXT \
			  VERSION.TXT

# Define how to build the documenation.

docdir			= $(top_builddir)/doc

apidoc:
	@echo Making $@ in $(docdir)
	@cd $(docdir) && make $@


# List the files that need to be cleaned.

MAINTAINERCLEANFILES	= Makefile.in \
			  aclocal.m4 \
			  configure \
			  config-h.in \
			  stamp-h.in \
			  $(AUX_DIST)
