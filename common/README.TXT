/**
@file README.TXT

@mainpage brick::common Utility Library

@section intro_sec Introduction

This library is the foundation for all of the other libraries in the
brick collection.  It provides classes for working with exceptions,
reference counting, portable numeric types with known precision, byte
order, and more.  The current version number is shown in file
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

You should be able to build this library without installing any
additional packages.

@section installation_sec Installation

@subsection autotools GNU Autotools

On every other platform, building and installing is done via GNU
Autotools.  If you're building from a tarball, the simplest
possible version looks like this:

@code
> tar -xvf brickcommon-x.y.z.tar.gz
> cd brickcommon-x.y.z
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
> [get source from cvs and run ./bootstrap, or unpack tarball]
> ./configure --prefix=/home/jfisher/software
> make
> make install
@endcode

Or some compiler options:

@code
> [get source from cvs and run ./bootstrap, or unpack tarball]
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
