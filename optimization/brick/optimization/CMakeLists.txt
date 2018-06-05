# Build file for the brickOptimization support library.

# This library is currently header-only, so no library target is needed.
#
# add_library(brickOptimization
#   foo.cc
#   )
# 
# install (TARGETS brickOptimization DESTINATION lib)

install (FILES

  autoGradientFunctionLM.hh
  gradientFunction.hh
  gradientFunctionLM.hh
  lossFunctions.hh
  optimizer.hh
  optimizerBFGS.hh
  optimizerCommon.hh
  optimizerLM.hh
  optimizerLineSearch.hh
  optimizerNelderMead.hh
  
  DESTINATION include/brick/optimization)

if (BRICK_BUILD_TESTS)
  add_subdirectory (test)
endif (BRICK_BUILD_TESTS)
