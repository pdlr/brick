/**
***************************************************************************
* @file iso12233/test/testImages.hh
*
* Header file specifying file names for images to be used in tests.
*
* Copyright (C) 2018 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_ISO12233_TEST_TESTIMAGES_HH
#define BRICK_ISO12233_TEST_TESTIMAGES_HH

#ifndef BRICK_TEST_DATA_DIR
#define BRICK_TEST_DATA_DIR "."
#endif

#include <string>
#include <brick/utilities/path.hh>

namespace brick {

  inline std::string
  getSlantedEdgeImageFileNamePNG0() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "slantedEdgeImage0.png");
  }

} // namespace brick

#endif /* #ifndef BRICK_ISO12233_TEST_TESTIMAGES_HH */
