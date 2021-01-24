# Build file for the brickIso12233 support library.

add_library(brickIso12233
  iso12233.cc
  )

# Instead of specifying -std=c++11 explicitly, we just tell CMake what
# features we need, and let it figure out the compiler flags.
target_compile_features(brickIso12233 PUBLIC
  cxx_long_long_type)

target_link_libraries (brickIso12233
  brickComputerVision
  brickLinearAlgebra
  brickNumeric
  )

install (TARGETS brickIso12233 DESTINATION lib)
install (FILES

  iso12233.hh iso12233_impl.hh
  oecf.hh oecf_impl.hh
  symbolImports.hh
  
  DESTINATION include/brick/iso12233)


add_subdirectory (demos)

if (BRICK_BUILD_TESTS)
  add_subdirectory (test)
endif (BRICK_BUILD_TESTS)
