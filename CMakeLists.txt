# Build file for the brick suite of support libraries.
#
# Recipes:
#  - To manually control installation directory:
#      cmake -DCMAKE_INSTALL_PREFIX=/foo/bar ..
#
#  - To prevent building tests:
#      cmake -DBRICK_BUILD_TESTS=OFF ..
#
#  - To build with debugging flags and no optimization:
#      cmake -DCMAKE_BUILD_TYPE=Debug ..
#
#  - To run tests:
#      make test
#
#  - To run tests and see output from any that fail:
#      env CTEST_OUTPUT_ON_FAILURE=1 make test
#
#  - To build the source package:
#      cmake -DCPACK_PACKAGE_VERSION=2.0.0 ..
#      make package_source
#
#  - To build brick as an external package in your own project, add this
#    (or something similar) to your CMakeLists.txt:
#      ExternalProject_Add(
#      brickLibraries
#      URL ${CMAKE_CURRENT_SOURCE_DIR}/external/Brick-2.0.0-Source.tar.gz
#      CMAKE_ARGS
#        -DBRICK_BUILD_SHARED_LIBRARIES=ON
#        -DBRICK_BUILD_TESTS=OFF
#        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_SOURCE_DIR}
#      )

cmake_minimum_required (VERSION 3.0)
project (Brick)

# Options for the project build.

option (BRICK_BUILD_SHARED_LIBRARIES "Controls whether generated libraries are static or shared." OFF)
option (BRICK_BUILD_TESTS "Build tests along with brick library code." ON)
option (BRICK_DEBUG_ARRAY_BOUNDS "Turn on run-time bounds checks." OFF)

option (BRICK_BUILD_COMMON
  "brickCommon is used by all other brick libraries."
  ON)
option (BRICK_BUILD_PORTABILITY
  "brickPortability provides platform independence to the brick library suite."
  ON)
option (BRICK_BUILD_TEST_LIB
  "brickTest is a simple unit testing framework.  You need this to test the other libraries in the brick library suite."
  ON)
option (BRICK_BUILD_NUMERIC
  "brickNumeric provides matrices, vectors, and basic math functions."
  ON)
option (BRICK_BUILD_PIXEL_GRAPHICS
  "brickPixelGraphics provides utilities for drawing directly on images."
  ON)
option (BRICK_BUILD_LINEAR_ALGEBRA
  "brickLinearAlgebra builds upon brickNumeric.  It needs lapack and blas."
  ON)
option (BRICK_BUILD_RANDOM
  "brickRandom provides simple pseudorandom number generators. It needs lapack and blas."
  ON)
option (BRICK_BUILD_OPTIMIZATION
  "brickOptimization provides simple nonlinear optimization routines."
  ON)
option (BRICK_BUILD_UTILITIES
  "brickUtilities contains things that didn't fit elsewhere, like argument parsing, string manipulation, etc."
  ON)
option (BRICK_BUILD_GEOMETRY
  "brickGeometry needs brickLinearAlgebra."
  ON)
option (BRICK_BUILD_COMPUTER_VISION
  "brickComputerVision depends on all of the other brick libraries."
  ON)
option (BRICK_BUILD_ISO12233
  "brickIso12233 depends on all of the other brick libraries."
  ON)

# Set up compiler flags.  Be as strict as possible.  The
# add_definitions line may not be portable.  We'll look into how to
# make it generic later.
# Note that we don't specify --std=c++11 here.  Instead, we use the 
# target_compile_features() function when setting the properties
# of the libraries themselves.

set (CMAKE_CXX_EXTENSIONS OFF)
set (BUILD_SHARED_LIBS ${BRICK_BUILD_SHARED_LIBRARIES})
add_definitions (
  -Wall
  -Wextra
  -Wshadow
  -Wpedantic
  -Werror)

if (BRICK_DEBUG_ARRAY_BOUNDS)
  add_definitions (-DBRICK_NUMERIC_CHECKBOUNDS=1)
endif (BRICK_DEBUG_ARRAY_BOUNDS)


# Other configuration.

if (BRICK_BUILD_TESTS)
  enable_testing()
endif (BRICK_BUILD_TESTS)

# Include each of the brick libraries.

set (BRICK_LIBS "")

macro (brick_configure_library directory_name library_name do_build)
  if (${do_build})
    include_directories ("${PROJECT_SOURCE_DIR}/${directory_name}")
    add_subdirectory (${directory_name})
    set (BRICK_LIBS ${BRICK_LIBS} ${library_name})
  endif (${do_build})
endmacro (brick_configure_library directory_name library_name do_build)

brick_configure_library (
  common brickCommon ${BRICK_BUILD_COMMON}
  )
brick_configure_library (
  portability brickPortability ${BRICK_BUILD_PORTABILITY}
  )
brick_configure_library (
  test brickTest ${BRICK_BUILD_TEST_LIB}
  )
brick_configure_library (
  numeric brickNumeric ${BRICK_BUILD_NUMERIC}
  )
brick_configure_library (
  pixelGraphics brickPixelGraphics ${BRICK_BUILD_PIXEL_GRAPHICS}
  )
brick_configure_library (
  linearAlgebra brickLinearAlgebra ${BRICK_BUILD_LINEAR_ALGEBRA}
  )
brick_configure_library (
  random brickRandom ${BRICK_BUILD_RANDOM}
  )
brick_configure_library (
  optimization brickOptimization ${BRICK_BUILD_OPTIMIZATION}
  )
brick_configure_library (
  utilities brickUtilities ${BRICK_BUILD_UTILITIES}
  )
brick_configure_library (
  geometry brickGeometry ${BRICK_BUILD_GEOMETRY}
  )
brick_configure_library (
  computerVision brickComputerVision ${BRICK_BUILD_COMPUTER_VISION}
  )
brick_configure_library (
  iso12233 brickIso12233 ${BRICK_BUILD_ISO12233}
  )


# Packaging.

# We'll need a way to package up releases, but first, a list of files
# we _don't_ want to include in the source packages.

# Don't package version control files.
set (CPACK_SOURCE_IGNORE_FILES
  ${CPACK_SOURCE_IGNORE_FILES}
  "^${CMAKE_CURRENT_SOURCE_DIR}/\\\\.hg.*"
  )

# Don't package build files.
set (CPACK_SOURCE_IGNORE_FILES
  ${CPACK_SOURCE_IGNORE_FILES}
  "^${CMAKE_CURRENT_SOURCE_DIR}/build.*"
  "scripts$"
  )

# Including CPack adds build targets that you can invoke with
# "make package" and "make package_source".
include (CPack)

# ========= Useful code that we don't need yet ============

# Here's how you link with brick libraries:
#
# # add_executable (foo foo.cc)
# # target_link_libraries (foo ${BRICK_LIBS})

# Useful code for setting version numbers.  Maybe we'll want this
# in subdirectories.
# 
# set (BRICK_NUMERIC_VERSION_MAJOR 1)
# set (BRICK_NUMERIC_VERSION_VERSION_MINOR 0)
# 
# # Propagate our version information to the source code.
# 
# configure_file (
#   "${PROJECT_SOURCE_DIR}/numeric/brick/version.hh.in"
#   "${PROJECT_BINARY_DIR}/numeric/brick/version.hh"
#   )
# 
# # Make sure the updated header file is visible via include
# # search path.
# 
# include_directories("${PROJECT_BINARY_DIR}/numeric")
#
# # Here's example content for version.hh.in:
# #
# #   #define BRICK_NUMERIC_VERSION_MAJOR @BRICK_NUMERIC_VERSION_MAJOR@
# #   #define BRICK_NUMERIC_VERSION_MINOR @BRICK_NUMERIC_VERSION_MINOR@

# Useful code for setting options.
#
# # option (MY_OPTION "Option description")
#
# And in version.hh.in:
#
# # #cmakedefine MY_OPTION
# # [...]
# # #ifdef MY_OPTION
# # [...]
# # #endif

