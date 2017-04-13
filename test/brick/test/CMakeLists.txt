# Build file for the brickTest support library.

add_library (brickTest autoregister.cc)
add_library (brickTestAutoMain automaticMain.cc)

install (TARGETS brickTest DESTINATION lib)
install (TARGETS brickTestAutoMain DESTINATION lib)
install (FILES

  autoregister.hh
  functors.hh
  runnableObject.hh
  testException.hh
  testFixture.hh
  testMacros.hh

  DESTINATION include/brick/test)


# This library doesn't have any tests yet.
# 
# if (BUILD_BRICK_TESTS)
#   add_subdirectory (test)
# endif (BUILD_BRICK_TESTS)
