# Build file for the brickCommon support library.

add_library(brickCommon
  compileTimestamp.cc
  exception.cc
  expect.cc
  traceable.cc
  )

# Instead of specifying -std=c++11 explicitly, we just tell CMake what
# features we need, and let it figure out the compiler flags.
target_compile_features(brickCommon PUBLIC
  cxx_long_long_type)

install (TARGETS brickCommon DESTINATION lib)
install (FILES

  argumentDescription.hh
  byteOrder.hh
  compileTimestamp.hh
  complexNumber.hh
  constants.hh
  exception.hh
  expect.hh
  functional.hh
  mathFunctions.hh
  referenceCount.hh
  stridedPointer.hh
  traceable.hh
  triple.hh
  types.hh
  versionNumber.hh
  
  DESTINATION include/brick/common)


if (BRICK_BUILD_TESTS)
  add_subdirectory (test)
endif (BRICK_BUILD_TESTS)
