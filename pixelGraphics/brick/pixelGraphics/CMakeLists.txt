# Build file for the brickPixelGraphics support library.

# This library is currently header-only, so no library target is needed.
# add_library(brickPixelGraphics
# 
#   )
# 
# install (TARGETS brickPixelGraphics DESTINATION bin)

install (FILES

  draw2D.hh draw2D_impl.hh
  
  DESTINATION include/brick/pixelGraphics)


# This library doesn't have any tests yet.
# 
# if (BRICK_BUILD_TESTS)
#   add_subdirectory (test)
# endif (BRICK_BUILD_TESTS)
