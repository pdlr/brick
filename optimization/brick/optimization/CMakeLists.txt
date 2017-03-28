# Build file for the brickOptimization support library.

# This library is currently header-only, so no library target is needed.
#
# add_library(brickOptimization
#   foo.cc
#   )
# 
# install (TARGETS brickOptimization DESTINATION bin)

install (FILES

  gradientFunction.hh
  gradientFunctionLM.hh
  optimizer.hh
  optimizerBFGS.hh
  optimizerCommon.hh
  optimizerLM.hh
  optimizerLineSearch.hh
  optimizerNelderMead.hh
  
  DESTINATION include/brick/optimization)


# This library doesn't have any tests yet.
# 
# if (BUILD_BRICK_TESTS)
#   add_subdirectory (test)
# endif (BUILD_BRICK_TESTS)
