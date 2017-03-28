# Build file for the brickCommon support library.

add_library(brickCommon
  compileTimestamp.cc
  exception.cc
  expect.cc
  traceable.cc
  )

install (TARGETS brickCommon DESTINATION bin)
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
