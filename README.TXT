/**
@file README.TXT

@mainpage Brick Utility Libraries

@section intro_sec Introduction

This is a collection of libraries I initially wrote when I was
consulting. I wanted a common platform of unencumbered code for use
with different clients, so I wouldn't have to start from scratch each
time. I use it as a foundation for much of my software.

The current version number is shown in file VERSION.TXT, and recent
changes can be seen in the file RELEASE_NOTES.TXT.  I've released it
under a pretty permissive license (LGPL) already, but let me know if
you want different terms.

@section platform_sec Platform

These libraries are in regular use under Linux (64-bit).  They should
be portable to any platform supporting C++11, but they haven't been
tested under windows in a very long time.

@section libraries_sec Libraries

There are several libraries in this collection, with the higher-level
libraries building on the lower-level ones.  Here is a quick ascii
table illustrating the dependencies between the libraries.  To save
columns, I've omitted the "brick" prefix in library names, and used
three-letter abbreviations in the column headers.  The final
rows/columns of this table capture dependencies on external libraries
(LAPACK/BLAS and libpng/libpng++).  We're currently using CMake to
build, and it should track these dependencies for you pretty
painlessly.

@verbatim
               Dependencies When Not Building Unit Tests
-----------------------------------------------------------------------------
       Library       |                   Depends on                          |
                     |CMN|PTB|TST|NUM|LAL|RND|OPT|UTL|GEO|CTV|PXG|MTF|LPK|PNG|  
Common         (CMN) | X |   |   |   |   |   |   |   |   |   |   |   |   |   |
Portability    (PTB) | X | X |   |   |   |   |   |   |   |   |   |   |   |   |
Test           (TST) | X |   | X |   |   |   |   |   |   |   |   |   |   |   |
Numeric        (NUM) | X |   |   | X |   |   |   |   |   |   |   |   |   |   |
LinearAlgebra  (LAL) | X |   |   | X | X |   |   |   |   |   |   |   | X |   |
Random         (RND) | X | X |   |   |   | X |   |   |   |   |   |   | X |   |
Optimization   (OPT) | X |   |   | X | X |   | X |   |   |   |   |   | X |   |
Utilities      (UTL) | X | X |   |   |   |   |   | X |   |   |   |   |   |   |
Geometry       (GEO) | X |   |   | X | X |   |   |   | X |   |   |   | X |   |
ComputerVision (CTV) | X |   |   | X | X | X | X |   | X | X |   |   | X | X |
PixelGraphics  (PXG) | X |   |   | X |   |   |   |   | X |   | X |   |   |   |
Iso12233       (MTF) | X |   |   | X | X | X | X |   | X | X |   | X | X | X |
LAPACK         (LPK) |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
libPng++       (PNG) |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
-----------------------------------------------------------------------------
@endverbatim

Each library comes with a bunch of unit tests.  If you want to build
the unit tests, the dependencies get a little deeper, as the tests
often depend on libraries that the tested code doesn't.  Again, CMake
should invisibly track this for you.

@verbatim
                  Dependencies When Building Unit Tests
-----------------------------------------------------------------------------
       Library       |                   Depends on                          |
                     |CMN|PTB|TST|NUM|LAL|RND|OPT|UTL|GEO|CTV|PXG|MTF|LPK|PNG|
Common         (CMN) | X |   |   |   |   |   |   |   |   |   |   |   |   |   |
Portability    (PTB) | X | X |   |   |   |   |   |   |   |   |   |   |   |   |
Test           (TST) | X |   | X |   |   |   |   |   |   |   |   |   |   |   |
Numeric        (NUM) | X | X | X | X |   |   |   |   |   |   |   |   |   |   |
LinearAlgebra  (LAL) | X |   | X | X | X |   |   |   |   |   |   |   | X |   |
Random         (RND) | X | X |   |   |   | X |   |   |   |   |   |   | X |   |
Optimization   (OPT) | X |   |   | X | X |   | X |   |   |   |   |   | X |   |
Utilities      (UTL) | X | X | X |   |   |   |   | X |   |   |   |   |   |   |
Geometry       (GEO) | X |   | X | X | X |   | X | X | X |   |   |   | X |   |
ComputerVision (CTV) | X | X | X | X | X | X | X | X | X | X |   |   | X | X |
PixelGraphics  (PXG) | X |   |   | X |   |   |   |   | X |   | X |   |   |   |
Iso12233       (MTF) | X | X | X | X | X | X | X | X | X | X |   | X | X | X |
LAPACK         (LPK) |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
libPng++       (PNG) |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
-----------------------------------------------------------------------------
@endverbatim

We group these libraries into three subgoups, as described in the following
subsections.

@subsection Core Libraries

The libraries listed in this section have been around for a
while. They are well tested, and their interfaces are stable. Features
get added from time to time, but we're generally careful not to break
existing user code.

@subsubsection brickCommon

This library is the foundation for all of the other libraries in the
brick collection.  It provides classes for working with exceptions,
reference counting (now largely obsoleted by shared pointers), portable
numeric types with known precision (the integer types are now obsoleted
by standard sized types, but the floating point types are still
relevant), byte order, and more.

@subsubsection brickPortability

This library was created to encapsulate whatever platform-specific
code is required by other, higher level libraries in the brick
suite.  Currently it contains some functions for dealing with
filesystem operations, some time-related functions, and an incomplete
implementation of snprintf.  For the most part, user code should not
access routines in brick::portability directly.  Most of this
functionality is already exposed through higher level libraries like
brickUtilities.  As with brickCommon, updates to the standard
have made a bunch of this code obsolete.

@subsubsection brickTest

This is a unit testing library, similar to CPPUNIT or boost::test, but
less featureful.  It provides classes and macros for writing tests and
coordinating their execution.  Currently doesn't include any facility
for interrupting tests that run too long, and lacks a formal schema
for test output, but it runs the brick tests without requiring any
outside dependencies.

@subsubsection brickNumeric

This is the original library in the collection. It has code for 1D,
2D, 3D, and ND arrays, coordinate transformations, and much more. This
is a good library to use if you need to port Numeric/numpy code from
Python to C++. Some recent additions to the library aren't stable, but
they're marked as such.

@subsubsection brickLinearAlgebra

This library provides C++ wrappers around LAPACK linear algebra
routines.  Currently it includes routines for matrix inverse, solution
of linear equations, SVD, determinant, eigenvectors/eigenvalues of
symmetric real matrices, etc.  A few routines are now templatized
to work with arbitrary types, but obviously you won't get the
LAPACK optimizations with types that aren't supported by LAPACK.

@subsubsection brickRandom

This library provides a C++ wrapper around the LAPACK pseudo-random
number generator.  It is is now largely obsoleted by updates to the
C++ standard library (see standard header <random>).

@subsubsection brickOptimization

This library provides classes implementing Levenberg-Marquardt
nonlinear least squares optimization, Quasi-Newton optimization, and
Downhill Simplex nonlinear optimization algorithms. These classes
follow the algorithms described in Numerical Recipes in C, by Press
and Flannery.  This has been a stable library for some time, but we
may extend it soon to take advantage of sparse matrices.  If this
happens, there may be some non-backwards compatible interface changes.
Note that the author now uses the Ceres nonlinear least squares
solver (ceres-solver.org) for most nonlinear least squares
optimization problems.

@subsubsection brickUtilities

This library provides some general routines & classes for manipulating
strings, parsing commandline arguments, doing filesystem operations,
working with dates & times, and a few other odds and ends.  It started
out as a C++ replacement for the python modules string and path.  It
is now gradually being replaced by updates to the C++ standard
library.

@subsection Unstable Libaries

The libraries listed in this section are still growing.  Please feel
free to use them, but bear in mind that they may change out from under you.

@subsubsection brickGeometry

This library provides some general routines for dealing with 2D and 3D
geometry. They originally lived in brickComputerVision (back when
brickComputerVision was called dlrComputerVision), but were broken out
into a separate library in hopes that they have value of their own.
It is still in the very early stages of development.

@subsubsection brickComputerVision

This library includes for basic morphological operations, edge
detection, simple deformable contours, histogram equalization, camera
models, and a continually growing list of computer vision related
things.  It's hard to compete with OpenCV these days, but the code
in this library is often more readable.

@subsubsection brickPixelGraphics

This library provides some classes and functions for drawing in images
and image-like objects.  It is still in early development.

@subsubsection brickIso12233

This library implements the e-SFR algorithm described in ISO 12233:2017.
Use it if you need to compute the MTF of an optical system using a
slanted edge target.  The code isn't super-robust yet (as of
May, 2019), but seems to be doing the right things.  It passes its
unit tests...

@subsubsection brickSparse

This library is intended to provide some classes and functions for
working with sparse matrices, but it hasn't been written yet.

@section contact_sec Contact

Bugfixes and patches welcome!  Please see the file LICENSE.TXT in
this directory for up-to-date contact information.

@section license_sec License

Please see the file LICENSE.TXT in this directory.

@section installation_sec Installation

The distribution includes CMakeLists.txt files to allow building with
CMake.  The traditional way to build on a linux machine is to start in
the top level of the source tree and issue these commands:

@code
> mkdir build
> cd build
> cmake -DCMAKE_INSTALL_PREFIX=/foo/bar ..
> make
> make test
> make install
@endcode

Thanks,
David LaRose
**/
