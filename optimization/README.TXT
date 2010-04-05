/**
@file README.TXT

@mainpage dlrOptimization Utility Library

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

The configure script for this package includes code generated from
Steven G. Johnson's autoconf files, acx_latex.m4 and acx_blas.m4.
Distribution of this generated code is permitted under the terms of
the associated license (GPLWithACException).  The actual files
acx_lapack.m4 and acx_blas.m4 should not be included in this
distribution.  If you would like access to these files, they are
available from http://ac-archive.sourceforge.net.

@section platform_sec Platform

This library is in regular use under Linux (32- & 64-bit), and
not-so-regular use under Windows (32-bit).  It should be portable
to any platform supporting ISO C++.

@section status_sec Status

This is a reasonably stable library.  Efforts will be made to maintain
a consistent interface, and not break existing code.

@section dependencies_sec Dependencies

This library currently consists only of header files.  To build using
these header files you will need to install install the following
libraries:

  dlrCommon
  dlrPortability
  dlrNumeric
  dlrLinearAlgebra
  lapack
  blas
  (and possibly, g2c)

In order to build the unit tests for this library, you must install
the dlrTest library.

You should be able to get dlrLinearAlgebra, dlrNumeric,
dlrPortability, and dlrCommon, from the same place you got this
library.  Lapack and blas are included with most linux distributions
and other unix-ish operating systems, such as OS X.  Most linux
distributions also include optimized versions of lapack and blas named
atlas or atlas3.  You can use these instead of the reference versions
for higher performance.  Non-free versions of lapack and blas which
are optimized for intel CPUs are available for various platforms from
http://www.intel.com.  You may be able to find pre-built lapack
binaries for your architecture at
http://www.netlib.org/lapack/archives/.  The library libg2c may be
needed if you compile with g++.  Libg2c is included as part of the gnu
compiler collection.  Depending on your OS distribution, it may be
prepackaged for you, in its own package, or bundled with g77.

@section installation_sec Installation

@subsection visual_studio Visual Studio

If you're building under Visual C++, load the solution file:

  ./visualc/dlrOptimization.sln

Building this file will put libraries in 

  "C:\Program Files\dlr_libs\{Release,Debug}"

and put include files in 

  "C:\Program Files\dlr_libs\include"

You can override these default installation locations from within
visual studio.  Note that if you do so, you may also have to update
the build files for other libraries in the dlr_libs suite so that they
know where to find the installed files from this library.

NOTE: Microsoft's compiler considers much of the C++ Standard Library
to be deprecated.  Building this code under recent versions of VC++
generates lots of warnings telling you to use the use the confusingly
named "Safe Standard C++ Library," which isn't standard at all, and as
far as I know isn't portable to any other platforms.  It's possible to
disable these warnings.  I encourage you to do so by adding
"/D_CRT_SECURE_NO_DEPRECATE" and "/D_SCL_SECURE_NO_DEPRECATE" to the
compiler command line via the "properties" dialog for the .vcproj file
within Visual Studio.  I've done this for the .vcproj files included
with dlr_libs, but you may find these warnings popping up in your own
code when you include dlr_libs header files.  Don't be fooled.

Please be aware that I hardly ever use windows, so the .sln/.vcproj
files are usually a little out of date.  I occasionally check that
friends who use these libraries under windows are still happy, but I
don't get patches too often.  Let me know if you have any trouble, and
please feel free to email with updates.  Also note that the .vcproj
files are all configured to build multi-threaded code.  This may cause
some trouble if you're linking with projects that use single threaded
code.  See the code generation page of the project properties dialog
to change this.

@subsection autotools GNU Autotools

On every other platform, building and installing is done via GNU
Autotools.  If you're building from a tarball, the simplest
possible version looks like this:

@code
> tar -xvf dlroptimization-x.y.z.tar.gz
> cd dlroptimization-x.y.z
> ./configure
> make
> make install
@endcode

If you're using the source files from subversion, you need to add a step:

@code
> [get source from svn]
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
