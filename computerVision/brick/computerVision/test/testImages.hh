/**
***************************************************************************
* @file computerVision/test/testImages.hh
*
* Header file specifying file names for images to be used in image
* processing tests.
*
* Copyright (C) 2005, 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_TEST_TESTIMAGES_HH
#define BRICK_COMPUTERVISION_TEST_TESTIMAGES_HH

#ifndef BRICK_TEST_DATA_DIR
#define BRICK_TEST_DATA_DIR "."
#endif

#include <string>
#include <brick/utilities/path.hh>

namespace brick {


  inline std::string
  getBullseyeFileNamePGM0() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "bullseyeTestImage0.pgm");
  }


  inline std::string
  getBullseyeFileNamePGM1() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "bullseyeTestImage1.pgm");
  }


  inline std::string
  getConnectedComponentsFileNamePGM0() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "connectedComponentsTestImage0.pgm");
  }


  inline std::string
  getConnectedComponentsFileNamePGM1() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "connectedComponentsTestImage1.pgm");
  }


  inline std::string
  getDilateErodeFileNamePGM0() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "dilateErodeTestImage0.pgm");
  }


  inline std::string
  getDilatedFileNamePGM0() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "dilatedTestImage0.pgm");
  }


  inline std::string
  getEdgeImageFileNamePGM0() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "edgeTestImage0.pgm");
  }


  inline std::string
  getErodedFileNamePGM0() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "erodedTestImage0.pgm");
  }


  inline std::string
  getTestImageFileNamePGM0() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "testImagePGM0.pgm");
  }


  inline std::string
  getTestImageFileNamePGM1() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "testImagePGM1.pgm");
  }


  inline std::string
  getTestImageFileNamePGM2() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "testImagePGM2.pgm");
  }


  inline std::string
  getTestImageFileNamePGM3() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "testImagePGM3.pgm");
  }


  inline std::string
  getTestImageFileNamePPM0() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "testImagePPM0.ppm");
  }


  inline std::string
  getTestImageFileNamePPM1() {
    return brick::utilities::joinPath(
      BRICK_TEST_DATA_DIR, "testImagePPM1.ppm");
  }

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_TEST_TESTIMAGES_HH */
