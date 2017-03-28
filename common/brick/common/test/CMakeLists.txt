include(CTest)

set (BRICK_COMMON_TEST_LIBS brickCommon)

# This macro simplifies building and adding test executables.

macro (brick_set_up_test test_name libraries)
  # Build the test in question.
  add_executable (${test_name} ${test_name}.cc)
  target_link_libraries (${test_name} ${libraries})

  # Arrange for the test to be run when the user executest the ctest command.
  add_test (${test_name}_target ${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_set_up_test test_name libraries)

# This macro specializes brick_set_up_test for brickCommon.

macro (brick_common_set_up_test test_name)
  brick_set_up_test (${test_name} ${BRICK_COMMON_TEST_LIBS})
endmacro (brick_common_set_up_test test_name)

# Here are the tests to be run.

brick_common_set_up_test (byteOrderTest)
brick_common_set_up_test (expectTest)
brick_common_set_up_test (referenceCountTest)
brick_common_set_up_test (traceableTest)
