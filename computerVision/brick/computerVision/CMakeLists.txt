# Build file for the brickComputerVision support library.

add_library(brickComputerVision
  connectedComponents.cc
  imageIO.cc
  histogramEqualize.cc
  keypointMatcherFast.cc
  keypointSelectorBullseye.cc
  keypointSelectorFast.cc
  pngReader.cc
  ransac.cc
  )

target_link_libraries (brickComputerVision
  brickLinearAlgebra
  brickNumeric
  ${PNG_LIBRARIES}
  )

if (PNG_FOUND)
  target_include_directories (brickComputerVision INTERFACE ${PNG_INCLUDE_DIRS})
endif (PNG_FOUND)

install (TARGETS brickComputerVision DESTINATION lib)
install (FILES

  calibrationTools.hh calibrationTools_impl.hh
  calibrationToolsRobust.hh calibrationToolsRobust_impl.hh
  cameraIntrinsics.hh
  cameraIntrinsicsDistortedPinhole.hh
  cameraIntrinsicsDistortedPinhole_impl.hh
  cameraIntrinsicsPinhole.hh cameraIntrinsicsPinhole_impl.hh
  cameraIntrinsicsPlumbBob.hh cameraIntrinsicsPlumbBob_impl.hh
  cameraIntrinsicsRational.hh cameraIntrinsicsRational_impl.hh
  canny.hh canny_impl.hh
  colorspaceConverter.hh colorspaceConverter_impl.hh
  connectedComponents.hh connectedComponents_impl.hh
  dilate.hh dilate_impl.hh
  disjointSet.hh disjointSet_impl.hh
  eightPointAlgorithm.hh eightPointAlgorithm_impl.hh
  erode.hh erode_impl.hh
  extendedKalmanFilter.hh extendedKalmanFilter_impl.hh
  featureAssociation.hh featureAssociation_impl.hh
  fitPolynomial.hh fitPolynomial_impl.hh
  fivePointAlgorithm.hh fivePointAlgorithm_impl.hh
  getEuclideanDistance.hh getEuclideanDistance_impl.hh
  histogramEqualize.hh
  image.hh
  imageFilter.hh imageFilter_impl.hh
  imageIO.hh imageIO_impl.hh
  imageFormat.hh
  imageFormatTraits.hh imageFormatTraits_impl.hh
  imagePyramid.hh imagePyramid_impl.hh
  imagePyramidBinomial.hh imagePyramidBinomial_impl.hh
  imageWarper.hh imageWarper_impl.hh
  kdTree.hh kdTree_impl.hh
  nTree.hh nTree_impl.hh
  kernel.hh kernel_impl.hh
  kernels.hh kernels_impl.hh
  keypointMatcherFast.hh keypointMatcherFast_impl.hh
  keypointSelectorBullseye.hh keypointSelectorBullseye_impl.hh
  keypointSelectorFast.hh keypointSelectorFast_impl.hh
  keypointSelectorHarris.hh keypointSelectorHarris_impl.hh
  naiveSnake.hh naiveSnake_impl.hh
  nonMaximumSuppress.hh nonMaximumSuppress_impl.hh
  nChooseKSampleSelector.hh nChooseKSampleSelector_impl.hh
  pixelBGRA.hh
  pixelHSV.hh
  pixelRGB.hh
  pixelRGBA.hh
  pixelYIQ.hh
  pngReader.hh
  randomSampleSelector.hh randomSampleSelector_impl.hh
  ransac.hh ransac_impl.hh
  ransacClassInterface.hh ransacClassInterface_impl.hh
  registerPoints3D.hh registerPoints3D_impl.hh
  segmenterFelzenszwalb.hh segmenterFelzenszwalb_impl.hh
  sobel.hh sobel_impl.hh
  stereoRectify.hh stereoRectify_impl.hh
  threePointAlgorithm.hh threePointAlgorithm_impl.hh
  thresholderSauvola.hh thresholderSauvola_impl.hh
  utilities.hh utilities_impl.hh
  
  DESTINATION include/brick/computerVision)


if (BRICK_BUILD_TESTS)
  add_subdirectory (test)
endif (BRICK_BUILD_TESTS)
