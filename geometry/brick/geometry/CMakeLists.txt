# Build file for the brickGeometry support library.

# This library is currently header-only, so no library target is needed.
#
# add_library(brickGeometry
#   foo.cc
#   )
# 
# install (TARGETS brickGeometry DESTINATION bin)

install (FILES

  bullseye2D.hh bullseye2D_impl.hh
  circle2D.hh circle2D_impl.hh
  circle3D.hh circle3D_impl.hh
  ellipse2D.hh ellipse2D_impl.hh
  lineSegment2D.hh lineSegment2D_impl.hh
  plane3D.hh plane3D_impl.hh
  ray2D.hh ray2D_impl.hh
  ray3D.hh ray3D_impl.hh
  triangle2D.hh triangle2D_impl.hh
  triangle3D.hh triangle3D_impl.hh
  utilities2D.hh utilities2D_impl.hh
  utilities3D.hh utilities3D_impl.hh

  DESTINATION include/brick/geometry)


if (BRICK_BUILD_TESTS)
  add_subdirectory (test)
endif (BRICK_BUILD_TESTS)
