/**
@file README.TXT

@mainpage brick::numeric Utility Library

@section intro_sec Introduction

This library provides some basic classes and functions for working
with numeric and scientific data.  It includes code for 1D, 2D, and 3D
arrays, coordinate transformations, and much more.  This is a good
library to use if you need to port Numpy code from Python to
C++.  The current version number is shown in file VERSION.TXT, and
recent changes can be seen in the file RELEASE_NOTES.TXT.

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

In order to build this library, you must first install the following
libraries:

  brick::common
  brick::portability

In order to build the unit tests for this library, you must install
the brick::test library.

You should be able to get these dependencies from the same place you
got this library.

@section installation_sec Installation

@subsection cmake CMake

The distribution now includes CMakeLists.txt files to allow building
with CMake.  For now, though, these files don't include a top-level
CMakeLists.txt.  The idea is that you might unpack the brick code into
your source tree and then include the brick CMakeLists.txt files from
your own CMake build structure.  If you want to build the brick
libraries stand-alone, there is an example top-level CMakeLists.txt
file available from the same place you got brickNumeric.

@subsection autotools GNU Autotools

Building and installing has traditionally been done via GNU Autotools.
If you're building from a tarball, the simplest possible version looks
like this:

@code
> tar -xvf bricknumeric-x.y.z.tar.gz
> cd bricknumeric-x.y.z
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
