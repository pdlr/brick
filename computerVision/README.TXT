/**
@file README.TXT

@mainpage brick::computerVision Utility Library

@section intro_sec Introduction

This library is still in the early stages of development.  Feel free
to play around with it, but please bear in mind that its interface is
not stable.  It includes code for basic morphological operations, edge
detection, simple deformable contours, histogram equalization, and
other computer vision related things.  The current version number is
shown in file VERSION.TXT, and recent changes can be seen in the file
RELEASE_NOTES.TXT.

@section contact_sec Contact

Bugfixes and patches welcome!  Please see the file LICENSE.TXT in
this directory for up-to-date contact information.

@section license_sec License

Please see the file LICENSE.TXT in this directory.

The configure script for this package includes code generated from
Steven G. Johnson's autoconf files, acx_latex.m4 and acx_blas.m4.
Distribution of this generated code is permitted under the terms of
the associated license (GPLWithACException).  The actual files
acx_lapack.m4 and acx_blas.m4 should not be included in this
distribution.  If you would like access to these files, they are
available from http://ac-archive.sourceforge.net.

@section platform_sec Platform

This library is in regular use under Linux (32- & 64-bit).  It should
be portable to any platform supporting ISO C++.

@section status_sec Status

This library is unstable.  Please expect its interface to change
significantly.

@section dependencies_sec Dependencies

In order to build this library, you must first install the following
libraries:

  brick::common
  brick::portability
  brick::numeric
  brick::linearAlgebra
  brick::geometry
  lapack
  blas
  (and possibly, g2c)

In order to build the unit tests for this library, you must install
the brick::test library.

You should be able to get brick::common, brick::portability, 
brick::numeric, brick::linearAlgebra, brick::geometry, and brick::test 
from the same place you got this library.  Lapack and blas are included 
with most linux distributions and other unix-ish operating systems, such 
as OS X.  Most linux distributions also include optimized versions of lapack 
and blas named atlas or atlas3.  You can use these instead of the reference
versions for higher performance.  Non-free versions of lapack and blas 
that are optimized for intel CPUs are available for various platforms 
from http://www.intel.com.  You may be able to find pre-built lapack
binaries for your architecture at  http://www.netlib.org/lapack/archives/.
The library libg2c may be needed if you compile with g++.Libg2c is included 
as part of the gnu compiler collection.  Depending on your OS distribution, 
it may be prepackaged for you, in its own package, or bundled with g77.

@section installation_sec Installation

@subsection autotools GNU Autotools

Building and installing is done via GNU Autotools.  If you're building 
from a tarball, the simplest possible version looks like this:

@code
> tar -xvf brickcomputervision-x.y.z.tar.gz
> cd brickcomputervision-x.y.z
> ./configure
> make
> make install
@endcode

If you're using the source files from mercurial, you need to add a step:

@code
> [get source from mercurial]
> ./bootstrap
> ./configure
> make
> make install
@endcode

Of course, you may want to specify an installation directory:

@code
> [get source from mercurial and run ./bootstrap, or unpack tarball]
> ./configure --prefix=/home/jfisher/software
> make
> make install
@endcode

Or some compiler options:

@code
> [get source from mercurial and run ./bootstrap, or unpack tarball]
> env CXXFLAGS="-g -Wall" ./configure \
    --prefix=/home/jfisher/software
> make
> make install
@endcode

For more information on configure options:

@code
> ./configure --help
@endcode

Thanks,
David LaRose
**/
