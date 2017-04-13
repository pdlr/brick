# Build file for the brickPortability support library.

add_library(brickPortability

  filesystem.cc
  standardC.cc
  timeUtilities.cc

  )

install (TARGETS brickPortability DESTINATION lib)
install (FILES

  filesystem.hh
  standardC.hh
  timeUtilities.hh
  
  DESTINATION include/brick/portability)


# This library doesn't have any tests yet.
#
# if (BUILD_BRICK_TESTS)
#   add_subdirectory (test)
# endif (BUILD_BRICK_TESTS)
