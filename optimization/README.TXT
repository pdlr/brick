/**
@file README.TXT

@mainpage brick::optimization Utility Library

@section intro_sec Introduction

This library provides classes implementing Levenberg-Marquardt
nonlinear least squares optimization, Quasi-Newton optimization, and
Downhill Simplex nonlinear optimization algorithms. These classes
follow the algorithms described in Numerical Recipes in C, by Press
and Flannery.  The current version number is shown in file
VERSION.TXT, and recent changes can be seen in the file
RELEASE_NOTES.TXT.

@section contact_sec Contact

Bugfixes and patches welcome!  Please see the file LICENSE.TXT in
this directory for up-to-date contact information.

@section license_sec License

Please see the file LICENSE.TXT in this directory.

@section platform_sec Platform

This library is in regular use under Linux (32- & 64-bit).  It should
be portable to any platform supporting ISO C++.

@section status_sec Status

This is a reasonably stable library.  Efforts will be made to maintain
a consistent interface, and not break existing code.

@section dependencies_sec Dependencies

This library currently consists only of header files.  To build using
these header files you will need to install install the following
libraries:

  brick::common
  brick::portability
  brick::numeric
  brick::linearAlgebra
  lapack
  blas
  (and possibly, g2c)

In order to build the unit tests for this library, you must install
the brick::test library.

You should be able to get brick::linearAlgebra, brick::numeric, 
brick::portability, brick::common, and brick::test from the same place you 
got this library.  Lapack and blas are included with most linux distributions 
and other unix-ish operating systems, such as OS X.  Most linux distributions 
also include optimized versions of lapack and blas named atlas or atlas3.  You 
can use these instead of the reference versions for higher performance.
Non-free versions of lapack and blas that are optimized for intel
CPUs are available for various platforms from http://www.intel.com.
You may be able to find pre-built lapack binaries for your
architecture at http://www.netlib.org/lapack/archives/.  The library
libg2c may be needed if you compile with g++.  Libg2c is included as
part of the gnu compiler collection.  Depending on your OS
distribution, it may be prepackaged for you, in its own package, or
bundled with g77.

@section installation_sec Installation

@subsection autotools GNU Autotools

@subsection cmake CMake

The distribution now includes CMakeLists.txt files to allow building
with CMake.  For now, though, these files don't include a top-level
CMakeLists.txt.  The idea is that you might unpack the brick code into
your source tree and then include the brick CMakeLists.txt files from
your own CMake build structure.  If you want to build the brick
libraries stand-alone, there is an example top-level CMakeLists.txt
file available from the same place you got brickOptimization.

@subsection autotools GNU Autotools

Building and installing has traditionally been done via GNU Autotools.
If you're building from a tarball, the simplest possible version looks
like this:

@code
> tar -xvf brickoptimization-x.y.z.tar.gz
> cd brickoptimization-x.y.z
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
