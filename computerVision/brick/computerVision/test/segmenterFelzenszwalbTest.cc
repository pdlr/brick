/**
***************************************************************************
* @file brick/computerVision/test/segmenterFelzenszwalbTest.cc
*
* Source file defining tests for routines declared in
* brick/computerVision/segmenterFelzenszwalb.hh.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/test/testImages.hh>
#include <brick/computerVision/segmenterFelzenszwalb.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/utilities.hh>
#include <brick/test/testFixture.hh>


using namespace brick::computerVision;

namespace brick {

  namespace computerVision {

    class SegmenterFelzenszwalbTest
      : public brick::test::TestFixture<SegmenterFelzenszwalbTest> {

    public:

      SegmenterFelzenszwalbTest();
      ~SegmenterFelzenszwalbTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testSegmenterFelzenszwalb();

    private:

    }; // class SegmenterFelzenszwalbTest


    /* ============== Member Function Definititions ============== */

    SegmenterFelzenszwalbTest::
    SegmenterFelzenszwalbTest()
      : brick::test::TestFixture<SegmenterFelzenszwalbTest>("SegmenterFelzenszwalbTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testSegmenterFelzenszwalb);
    }


    void
    SegmenterFelzenszwalbTest::
    testSegmenterFelzenszwalb()
    {
      SegmenterFelzenszwalb<EdgeDefaultFunctor<double>, double> segmenter(200, 0.8, 20);
      Image<GRAY8> inputImage0 = readPGM8(getTestImageFileNamePGM0());
      segmenter.segment(inputImage0);
      brick::numeric::Array2D<brick::common::UnsignedInt32> labelArray =
        segmenter.getLabelArray();
      brick::numeric::Array2D<brick::common::UnsignedInt16> labelImage(
        labelArray.rows(), labelArray.columns());
      labelImage.copy(labelArray);
      writePGM16("foo.pgm", labelImage);
    }

  } // namespace computerVision

} // namespace brick

#if 0

int main(int argc, char** argv)
{
  brick::computerVision::SegmenterFelzenszwalbTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::SegmenterFelzenszwalbTest currentTest;

}

#endif
