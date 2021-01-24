include(CTest)

set (BRICK_ISO12233_DEMO_LIBS
  brickIso12233
  brickComputerVision
  brickUtilities
  brickPortability
  brickCommon
  )


# This macro simplifies building and adding demo executables.

macro (brick_iso12233_set_up_demo demo_name)
  # Build the demo in question.
  add_executable (iso12233_${demo_name} ${demo_name}.cc)
  target_link_libraries (iso12233_${demo_name} ${BRICK_ISO12233_DEMO_LIBS})
endmacro (brick_iso12233_set_up_demo demo_name)


# Here are the demos to be built.

brick_iso12233_set_up_demo (analyzePatch)
