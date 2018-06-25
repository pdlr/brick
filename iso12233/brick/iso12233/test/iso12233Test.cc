/**
***************************************************************************
* @file brick/iso12233/test/iso12233Test.cc
*
* Source file defining tests for the MTF calculating routines.
*
* Copyright (C) 2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/iso12233/iso12233.hh>

#include <brick/test/testFixture.hh>
#include <brick/utilities/timeUtilities.hh>

namespace brick {

  namespace iso12233 {

    class Iso12233Test
      : public brick::test::TestFixture<Iso12233Test> {

    public:

      Iso12233Test();
      ~Iso12233Test() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testIso12233();

    private:

      double m_defaultTolerance;

    }; // class Iso12233Test


    /* ============== Member Function Definititions ============== */

    Iso12233Test::
    Iso12233Test()
      : brick::test::TestFixture<Iso12233Test>("Iso12233Test"),
        m_defaultTolerance(1.0E-8)
    {
      BRICK_TEST_REGISTER_MEMBER(testIso12233);
    }


    void
    Iso12233Test::
    testIso12233()
    {
      constexpr std::size_t patchWidth = 100;
      constexpr std::size_t patchHeight = 100;
      constexpr std::size_t windowWidth = 50;

      // Create a test image with a vertically straight up and down
      // edge.  We should have trouble with this image because the
      // vertical edge doesn't let us so super-resolution.
      Image<brick::computerVision::GRAY8> edgeImage(
        patchHeight, patchWidth); // Rows, columns.
      for(std::size_t rr = 0; rr < patchHeight; ++rr) {
        std::size_t cc = 0;
        while(cc < ((patchWidth / 2) - (windowWidth / 4))) {
          edgeImage(rr, cc) = 100;
          ++cc;
        }
        while(cc < patchWidth) {
          edgeImage(rr, cc) = 200;
          ++cc;
        }
      }

      // Try to process the image.
      Array1D<double> mtf;
      BRICK_TEST_ASSERT_EXCEPTION(
        brick::common::ValueException,
        mtf = iso12233<double>(edgeImage, windowWidth,
                               [](double arg){return arg;}));
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(/* int argc, char** argv */)
{
  brick::computerVision::Iso12233Test currentTest;
  bool result = currentTest.run();

  currentTest.exerciseKeypointSelectorFast("testImagePGM0.pgm");

  return (result ? 0 : 1);
}

#else

namespace {

  brick::iso12233::Iso12233Test currentTest;

}

#endif
