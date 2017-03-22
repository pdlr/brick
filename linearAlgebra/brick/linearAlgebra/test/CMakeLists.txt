include(CTest)

set (BRICK_LINEAR_ALGEBRA_TEST_LIBS

  brickLinearAlgebra
  brickNumeric
  brickPortability
  brickCommon
  brickTest
  brickTestAutoMain
  lapack
  blas

  )

# This macro simplifies building and adding test executables.

macro (brick_linear_algebra_set_up_test test_name)
  # Build the test in question.
  add_executable (${test_name} ${test_name}.cc)
  target_link_libraries (${test_name} ${BRICK_LINEAR_ALGEBRA_TEST_LIBS})

  # Arrange for the test to be run when the user executest the ctest command.
  add_test (${test_name}_target ${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_linear_algebra_set_up_test test_name)

# Here are all the tests to be run.

brick_linear_algebra_set_up_test(linearAlgebraTest)
