# Build file for the brickRandom support library.

add_library(brickRandom

  pseudoRandom.cc

  )

install (TARGETS brickRandom DESTINATION bin)
install (FILES

  clapack.hh
  pseudoRandom.hh
  
  DESTINATION include/brick/random)

# libBrickRandom doesn't have any tests yet.
# 
# if (BRICK_BUILD_TESTS)
#   add_subdirectory (test)
# endif (BRICK_BUILD_TESTS)
