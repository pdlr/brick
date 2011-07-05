#
# Makefile.am
#   Process this file with automake to produce Makefile.in.
#
# $Revision: 842 $
# $Date: 2006-12-15 01:09:48 -0500 (Fri, 15 Dec 2006) $
#


# Set VPATH so we can build in a different directory.

srcdir			= @srcdir@
VPATH			= @srcdir@


#

doxygen_config		= $(srcdir)/doxygen.config

apidoc:
	doxygen $(doxygen_config)
	cd latex && sed "s/\\setlength{\\footrulewidth}{/\\renewcommand{\\footrulewidth}{/g" < refman.tex > refman2.tex && mv refman2.tex refman.tex

# Specify the list of additional files that need to be included in the
# distribution.

EXTRA_DIST		= $(doxygen_config)


# List the files that need to be cleaned.

CLEANFILES		= $(wildcard html/*) $(wildcard latex/*)

MAINTAINERCLEANFILES	= Makefile.in

