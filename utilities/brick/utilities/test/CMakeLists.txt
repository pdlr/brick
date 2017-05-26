include(CTest)

set (BRICK_UTILITIES_TEST_LIBS
  brickUtilities
  brickTest
  brickTestAutoMain
  )

# This macro simplifies building and adding test executables.

macro (brick_utilities_set_up_test test_name)
  # Build the test in question.
  add_executable (utilities_${test_name} ${test_name}.cc)
  target_link_libraries (utilities_${test_name} ${BRICK_UTILITIES_TEST_LIBS})

  # Arrange for the test to be run when the user executest the ctest command.
  add_test (utilities_${test_name}_target utilities_${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (utilities_${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_utilities_set_up_test test_name)

# Here are the tests to be run.

brick_utilities_set_up_test (dateTest)
brick_utilities_set_up_test (imageIOTest)
brick_utilities_set_up_test (lockFileTest)
brick_utilities_set_up_test (optionParserTest)
brick_utilities_set_up_test (pathTest)
brick_utilities_set_up_test (stringManipulationTest)
brick_utilities_set_up_test (teeTest)
brick_utilities_set_up_test (timeUtilitiesTest)
