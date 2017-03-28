# Build file for the brickNumeric support library.

add_library(brickNumeric

  ieeeFloat32.cc
  index2D.cc
  index3D.cc
  transform2D.cc
  transform3D.cc
  transform3DTo2D.cc

  )

target_link_libraries (brickNumeric
  brickCommon
  )

install (TARGETS brickNumeric DESTINATION bin)
install (FILES

  amanatidesWoo2D.hh amanatidesWoo2D_impl.hh
  amanatidesWoo2DIterator.hh amanatidesWoo2DIterator_impl.hh
  amanatidesWoo3D.hh amanatidesWoo3D_impl.hh
  amanatidesWoo3DIterator.hh amanatidesWoo3DIterator_impl.hh
  array1D.hh array1D_impl.hh
  array2D.hh array2D_impl.hh
  array3D.hh array3D_impl.hh
  arrayND.hh arrayND_impl.hh
  bilinearInterpolator.hh bilinearInterpolator_impl.hh
  boxIntegrator2D.hh boxIntegrator2D_impl.hh
  bSpline.hh bSpline_impl.hh
  bSpline2D.hh bSpline2D_impl.hh
  convolutionStrategy.hh
  convolve1D.hh convolve1D_impl.hh
  convolve2D.hh convolve2D_impl.hh
  convolveND.hh convolveND_impl.hh
  derivativeRidders.hh derivativeRidders_impl.hh
  differentiableScalar.hh differentiableScalar_impl.hh
  extremumRecorder.hh
  functional.hh
  ieeeFloat32.hh
  index2D.hh
  index3D.hh
  filter.hh filter_impl.hh
  geometry2D.hh geometry2D_impl.hh
  mathFunctions.hh
  maxRecorder.hh
  minRecorder.hh
  normalizedCorrelator.hh normalizedCorrelator_impl.hh
  numericTraits.hh
  polynomial.hh polynomial_impl.hh
  quaternion.hh quaternion_impl.hh
  rotations.hh rotations_impl.hh
  slice.hh
  sampledFunctions.hh sampledFunctions_impl.hh
  scatteredDataInterpolator2D.hh scatteredDataInterpolator2D_impl.hh
  solveCubic.hh solveCubic_impl.hh
  solveQuadratic.hh solveQuadratic_impl.hh
  solveQuartic.hh solveQuartic_impl.hh
  staticArray1D.hh staticArray1D_impl.hh
  stencil2D.hh stencil2D_impl.hh
  subArray1D.hh subArray1D_impl.hh
  subArray2D.hh subArray2D_impl.hh
  subpixelInterpolate.hh subpixelInterpolate_impl.hh
  transform2D.hh transform2D_impl.hh
  transform3D.hh transform3D_impl.hh
  transform3DTo2D.hh transform3DTo2D_impl.hh
  typeConverter.hh
  vector2D.hh vector2D_impl.hh
  vector3D.hh vector3D_impl.hh
  utilities.hh utilities_impl.hh
  
  DESTINATION include/brick/numeric)


if (BRICK_BUILD_TESTS)
  add_subdirectory (test)
endif (BRICK_BUILD_TESTS)
