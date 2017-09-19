include(CTest)

set (BRICK_ISO12233_TEST_LIBS
  brickIso12233
  brickComputerVision
  brickUtilities
  brickRandom
  brickPortability
  brickTest
  brickTestAutoMain
  )

# This macro simplifies building and adding test executables.

macro (brick_iso12233_set_up_test test_name)
  # Build the test in question.
  add_executable (iso12233_${test_name} ${test_name}.cc)
  target_link_libraries (iso12233_${test_name} ${BRICK_ISO12233_TEST_LIBS})

  # Arrange for the test to be run when the user executes the ctest command.
  add_test (iso12233_${test_name}_target iso12233_${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_iso12233_set_up_test test_name)

# This macro simplifies copying support files to the directory where
# tests will be run.

macro (brick_iso12233_set_up_test_data file_name)
  configure_file (
    ${file_name}
    ${CMAKE_CURRENT_BINARY_DIR}/${file_name}
    COPYONLY)
endmacro (brick_iso12233_set_up_test_data file_name)


# Here are the tests to be run.

brick_iso12233_set_up_test (iso12233Test)
