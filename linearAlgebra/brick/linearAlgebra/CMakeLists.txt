# Build file for the brickLinearAlgebra support library.

# For when we get slick enough to track version information inside
# our CMakeLists.txt files.
# 
# set (BRICK_LINEARALGEBRA_VERSION 1.2)

add_library(brickLinearAlgebra
  linearAlgebra.cc
  )

target_link_libraries (brickLinearAlgebra
  brickCommon
  brickNumeric
  ${LAPACK_LIBRARIES}
  )

install (TARGETS brickLinearAlgebra DESTINATION lib)
install (FILES
  clapack.hh
  linearAlgebra.hh linearAlgebra_impl.hh
  DESTINATION include/brick/linearAlgebra)


# For when we get slick enough to track version information inside
# our CMakeLists.txt files.
# 
# set_target_properties (brickLinearAlgebra PROPERTIES
#    VERSION ${BRICK_LINEARALGEBRA_VERSION}
#    SOVERSION ${BRICK_LINEARALGEBRA_VERSION}
# )
  
if (BRICK_BUILD_TESTS)
  add_subdirectory (test)
endif (BRICK_BUILD_TESTS)
