/**
***************************************************************************
* @file brick/utilities/test/imageIOTest.cc
*
* Source file defining tests for image input and output routines.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <sstream>
#include <string>

#include <brick/utilities/imageIO.hh>
#include <brick/test/testFixture.hh>

using brick::test::TestFixture;

namespace brick {

  namespace utilities {
    
    class ImageIOTest : public TestFixture<ImageIOTest> {

    public:

      ImageIOTest();
      ~ImageIOTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Tests of non-member functions.
      void testWritePGM();
      void testWritePPM();

    }; // class ImageIOTest


    /* ============== Member Function Definititions ============== */

    ImageIOTest::
    ImageIOTest()
      : TestFixture<ImageIOTest>("ImageIOTest")
    {
      // Tests of non-member functions.
      BRICK_TEST_REGISTER_MEMBER(testWritePGM);
      BRICK_TEST_REGISTER_MEMBER(testWritePPM);
    }


    void
    ImageIOTest::
    testWritePGM()
    {
      const int imageRows = 480;
      const int imageColumns = 640;

      unsigned char* imageBuffer0 =
        new unsigned char[imageRows * imageColumns];
      imageBuffer0[0] = 1;
      writePGM("foo.pgm", imageBuffer0, imageRows, imageColumns,
               true, false, 8);
      delete[] imageBuffer0;

      unsigned short* imageBuffer1 =
        new unsigned short[imageRows * imageColumns];
      imageBuffer1[0] = 1;
      writePGM("bar.pgm", imageBuffer1, imageRows, imageColumns,
               true, false, 16);
      delete[] imageBuffer1;

      // Warning(xxx): Incomplete test.

      std::remove("foo.pgm");
      std::remove("bar.pgm");
    }


    void
    ImageIOTest::
    testWritePPM()
    {
      const int imageRows = 480;
      const int imageColumns = 640;

      unsigned char* imageBuffer0 =
        new unsigned char[imageRows * imageColumns * 3];
      imageBuffer0[0] = 1;
      writePPM("foo.ppm", imageBuffer0, imageRows, imageColumns,
               true, false, 8);
      delete[] imageBuffer0;

      unsigned short* imageBuffer1 =
        new unsigned short[imageRows * imageColumns * 3];
      imageBuffer1[0] = 1;
      writePPM("bar.ppm", imageBuffer1, imageRows, imageColumns,
               true, false, 16);
      delete[] imageBuffer1;

      // Warning(xxx): Incomplete test.

      std::remove("foo.ppm");
      std::remove("bar.ppm");
    }

  } // namespace utilities
  
} // namespace brick



#ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION

int main(int argc, char** argv)
{
  brick::utilities::ImageIOTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */

namespace {

  brick::utilities::ImageIOTest currentTest;
  
}

#endif /* #ifdef BRICK_TEST_NO_AUTOMATIC_REGISTRATION */
