# Build file for the brickComputerVision support library.

find_package (PNG)
if (PNG_FOUND)
    find_path (PNGPP_INCLUDE_DIR png++/png.hpp)
    if (NOT PNGPP_INCLUDE_DIR)
      add_definitions("-DHAVE_LIBPNG=0")
      set(PNGPP_INCLUDE_DIR "")
    endif ()
    
    set (PNG_INCLUDE_DIRS ${PNG_INCLUDE_DIRS} ${PNGPP_INCLUDE_DIR})
else (PNG_FOUND)
  add_definitions("-DHAVE_LIBPNG=0")
  set(PNGPP_INCLUDE_DIR "")
endif ()

add_subdirectory (brick/computerVision)
