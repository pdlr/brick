# Build file for the brickUtilities support library.

add_library(brickUtilities
  
  date.cc
  lockFile.cc
  optionDescription.cc
  optionParser.cc
  path.cc
  stringManipulation.cc

  )

target_link_libraries (brickUtilities
  brickPortability
  brickCommon
  )

install (TARGETS brickUtilities DESTINATION lib)
install (FILES

  date.hh
  exception.hh
  imageIO.hh
  frequencyGoverner.hh
  lockFile.hh
  optionDescription.hh
  optionParser.hh
  path.hh
  pythonIO.hh
  stringManipulation.hh
  timeUtilities.hh

  DESTINATION include/brick/utilities)


if (BRICK_BUILD_TESTS)
  add_subdirectory (test)
endif (BRICK_BUILD_TESTS)
