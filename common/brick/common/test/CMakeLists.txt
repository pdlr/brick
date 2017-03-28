include(CTest)

set (BRICK_COMMON_TEST_LIBS brickCommon)

# This macro simplifies building and adding test executables.

macro (brick_common_set_up_test test_name)
  # Build the test in question.
  add_executable (common_${test_name} ${test_name}.cc)
  target_link_libraries (common_${test_name} ${BRICK_COMMON_TEST_LIBS})

  # Arrange for the test to be run when the user executest the ctest command.
  add_test (common_${test_name}_target common_${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (common_${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_common_set_up_test test_name libraries)

# Here are the tests to be run.

brick_common_set_up_test (byteOrderTest)
brick_common_set_up_test (expectTest)
brick_common_set_up_test (referenceCountTest)
brick_common_set_up_test (traceableTest)
