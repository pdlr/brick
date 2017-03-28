include(CTest)

set (BRICK_COMPUTER_VISION_TEST_LIBS
  brickComputerVision
  brickUtilities
  brickRandom
  brickPortability
  brickTest
  brickTestAutoMain
  )

# This macro simplifies building and adding test executables.

macro (brick_computer_vision_set_up_test test_name)
  # Build the test in question.
  add_executable (computerVision_${test_name} ${test_name}.cc)
  target_link_libraries (computerVision_${test_name} ${BRICK_COMPUTER_VISION_TEST_LIBS})

  # Arrange for the test to be run when the user executest the ctest command.
  add_test (computerVision_${test_name}_target computerVision_${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_computer_vision_set_up_test test_name)

# Here are the tests to be run.

brick_computer_vision_set_up_test (calibrationToolsTest)
brick_computer_vision_set_up_test (calibrationToolsRobustTest)
brick_computer_vision_set_up_test (cameraIntrinsicsPinholeTest)
brick_computer_vision_set_up_test (cameraIntrinsicsPlumbBobTest)
brick_computer_vision_set_up_test (cameraIntrinsicsRationalTest)
brick_computer_vision_set_up_test (cannyTest)
brick_computer_vision_set_up_test (colorspaceConverterTest)
brick_computer_vision_set_up_test (connectedComponentsTest)
brick_computer_vision_set_up_test (dilateTest)
brick_computer_vision_set_up_test (eightPointAlgorithmTest)
brick_computer_vision_set_up_test (erodeTest)
brick_computer_vision_set_up_test (extendedKalmanFilterTest)
brick_computer_vision_set_up_test (featureAssociationTest)
brick_computer_vision_set_up_test (fivePointAlgorithmTest)
brick_computer_vision_set_up_test (getEuclideanDistanceTest)
brick_computer_vision_set_up_test (imageFilterTest)
brick_computer_vision_set_up_test (imageIOTest)
brick_computer_vision_set_up_test (imagePyramidTest)
brick_computer_vision_set_up_test (imagePyramidBinomialTest)
brick_computer_vision_set_up_test (imageWarperTest)
brick_computer_vision_set_up_test (kdTreeTest)
brick_computer_vision_set_up_test (keypointMatcherFastTest)
brick_computer_vision_set_up_test (keypointSelectorBullseyeTest)
brick_computer_vision_set_up_test (keypointSelectorFastTest)
brick_computer_vision_set_up_test (keypointSelectorHarrisTest)
brick_computer_vision_set_up_test (naiveSnakeTest)
brick_computer_vision_set_up_test (nChooseKSampleSelectorTest)
brick_computer_vision_set_up_test (nonMaximumSuppressTest)
# brick_computer_vision_set_up_test (ransacTest)
brick_computer_vision_set_up_test (registerPoints3DTest)
brick_computer_vision_set_up_test (segmenterFelzenszwalbTest)
brick_computer_vision_set_up_test (sobelTest)
brick_computer_vision_set_up_test (stereoRectifyTest)
brick_computer_vision_set_up_test (threePointAlgorithmTest)
brick_computer_vision_set_up_test (thresholderSauvolaTest)
brick_computer_vision_set_up_test (utilitiesTest)
