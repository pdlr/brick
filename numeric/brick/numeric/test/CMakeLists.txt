include(CTest)

set (BRICK_NUMERIC_TEST_LIBS
  brickNumeric
  brickPortability
  brickTest
  brickTestAutoMain
  )

# This macro simplifies building and adding test executables.

macro (brick_numeric_set_up_test test_name)
  # Build the test in question.
  add_executable (${test_name} ${test_name}.cc)
  target_link_libraries (${test_name} ${BRICK_NUMERIC_TEST_LIBS})

  # Arrange for the test to be run when the user executest the ctest command.
  add_test (${test_name}_target ${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_numeric_set_up_test test_name)

# Here are the tests to be run.

brick_numeric_set_up_test(amanatidesWoo2DTest)
brick_numeric_set_up_test(amanatidesWoo3DTest)
brick_numeric_set_up_test(array1DTest)
brick_numeric_set_up_test(array2DTest)
brick_numeric_set_up_test(array3DTest)
brick_numeric_set_up_test(arrayNDTest)
brick_numeric_set_up_test(bilinearInterpolatorTest)
brick_numeric_set_up_test(boxIntegrator2DTest)
brick_numeric_set_up_test(bSplineTest)
brick_numeric_set_up_test(bSpline2DTest)
brick_numeric_set_up_test(convolve1DTest)
brick_numeric_set_up_test(convolve2DTest)
brick_numeric_set_up_test(convolveNDTest)
brick_numeric_set_up_test(derivativeRiddersTest)
brick_numeric_set_up_test(differentiableScalarTest)
brick_numeric_set_up_test(filterTest)
brick_numeric_set_up_test(ieeeFloat32Test)
brick_numeric_set_up_test(index3DTest)
brick_numeric_set_up_test(geometry2DTest)
brick_numeric_set_up_test(maxRecorderTest)
brick_numeric_set_up_test(minRecorderTest)
brick_numeric_set_up_test(normalizedCorrelatorTest)
brick_numeric_set_up_test(polynomialTest)
brick_numeric_set_up_test(rotationsTest)
brick_numeric_set_up_test(sampledFunctionsTest)
brick_numeric_set_up_test(scatteredDataInterpolator2DTest)
brick_numeric_set_up_test(solveCubicTest)
brick_numeric_set_up_test(solveQuadraticTest)
brick_numeric_set_up_test(solveQuarticTest)
brick_numeric_set_up_test(stencil2DTest)
brick_numeric_set_up_test(subpixelInterpolateTest)
brick_numeric_set_up_test(transform2DTest)
brick_numeric_set_up_test(transform3DTest)
brick_numeric_set_up_test(utilitiesTest)
