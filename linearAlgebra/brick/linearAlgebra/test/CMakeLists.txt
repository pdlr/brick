include(CTest)

set (BRICK_LINEAR_ALGEBRA_TEST_LIBS
  brickLinearAlgebra
  brickTest
  brickTestAutoMain
  )

# This macro simplifies building and adding test executables.

macro (brick_linear_algebra_set_up_test test_name)
  # Build the test in question.
  add_executable (linearAlgebra_${test_name} ${test_name}.cc)
  target_link_libraries (linearAlgebra_${test_name} ${BRICK_LINEAR_ALGEBRA_TEST_LIBS})

  # Arrange for the test to be run when the user executest the ctest command.
  add_test (linearAlgebra_${test_name}_target linearAlgebra_${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (linearAlgebra_${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_linear_algebra_set_up_test test_name)

# Here are all the tests to be run.

brick_linear_algebra_set_up_test(linearAlgebraTest)
