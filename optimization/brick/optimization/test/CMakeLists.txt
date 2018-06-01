include(CTest)

set (BRICK_OPTIMIZATION_TEST_LIBS
  brickLinearAlgebra
  brickTest
  brickTestAutoMain
  )

# This macro simplifies building and adding test executables.

macro (brick_optimization_set_up_test test_name)
  # Build the test in question.
  add_executable (optimization_${test_name} ${test_name}.cc)
  target_link_libraries (optimization_${test_name} ${BRICK_OPTIMIZATION_TEST_LIBS})

  # Arrange for the test to be run when the user executest the ctest command.
  add_test (optimization_${test_name}_target optimization_${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (optimization_${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_optimization_set_up_test test_name)

# Here are all the tests to be run.

brick_optimization_set_up_test (autoGradientFunctionLMTest)
brick_optimization_set_up_test (lossFunctionsTest)
