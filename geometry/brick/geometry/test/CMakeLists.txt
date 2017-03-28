include(CTest)

set (BRICK_GEOMETRY_TEST_LIBS
  brickLinearAlgebra
  brickTest
  brickTestAutoMain
  )

# This macro simplifies building and adding test executables.

macro (brick_geometry_set_up_test test_name)
  # Build the test in question.
  add_executable (geometry_${test_name} ${test_name}.cc)
  target_link_libraries (geometry_${test_name} ${BRICK_GEOMETRY_TEST_LIBS})

  # Arrange for the test to be run when the user executest the ctest command.
  add_test (geometry_${test_name}_target geometry_${test_name})

  # All brick unit tests return 0 on success, nonzero otherwise,
  # so no need to set special properties that catch failures.
  # 
  # # set_tests_properties (geometry_${test_name}_target
  # #   PROPERTIES PASS_REGULAR_EXPRESSION "All tests pass")
endmacro (brick_geometry_set_up_test test_name)

# Here are all the tests to be run.

brick_geometry_set_up_test (bullseye2DTest)
brick_geometry_set_up_test (circle3DTest)
brick_geometry_set_up_test (ellipse2DTest)
brick_geometry_set_up_test (plane3DTest)
brick_geometry_set_up_test (ray2DTest)
brick_geometry_set_up_test (triangle2DTest)
brick_geometry_set_up_test (triangle3DTest)
brick_geometry_set_up_test (utilities2DTest)
brick_geometry_set_up_test (utilities3DTest)
