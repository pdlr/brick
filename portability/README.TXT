/**
@file README.TXT

@mainpage brick::portability Utility Library

@section intro_sec Introduction

This library was created to encapsulate whatever platform-specific
code is required by other, higher level libraries in the brick
suite.  Currently it contains some functions for dealing with
filesystem operations, some time-related functions, and an incomplete
implementation of snprintf.  For the most part, user code should not
access routines in brick::portability directly.  Most of this
functionality is already exposed through higher level libraries like
brick::utilities.  The current version number is shown in file
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

In order to build brick::portability, you must first install the brick::common
library.  You should be able to get brick::common from the same place you
got brick::portability.

@section installation_sec Installation

@subsection autotools GNU Autotools

Building and installing is done via GNU Autotools.  If you're building 
from a tarball, the simplest possible version looks like this:

@code
> tar -xvf brickportability-x.y.z.tar.gz
> cd brickportability-x.y.z
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
