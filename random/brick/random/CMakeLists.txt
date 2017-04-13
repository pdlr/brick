# Build file for the brickRandom support library.

add_library(brickRandom
  pseudoRandom.cc
  )

# Instead of specifying -std=c++11 explicitly, we just tell CMake what
# features we need, and let it figure out the compiler flags.
target_compile_features(brickRandom PUBLIC
  cxx_long_long_type)

install (TARGETS brickRandom DESTINATION lib)
install (FILES
  clapack.hh
  pseudoRandom.hh
  DESTINATION include/brick/random)

# libBrickRandom doesn't have any tests yet.
# 
# if (BRICK_BUILD_TESTS)
#   add_subdirectory (test)
# endif (BRICK_BUILD_TESTS)
