# Build file for the brickLinearAlgebra support library.

add_library(brickLinearAlgebra

  linearAlgebra.cc
  
  )

install (TARGETS brickLinearAlgebra DESTINATION bin)
install (FILES

  clapack.hh
  linearAlgebra.hh
  
  DESTINATION include)


if (BUILD_BRICK_TESTS)
  add_subdirectory (test)
endif (BUILD_BRICK_TESTS)
