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

#include <brick/numeric/bilinearInterpolator.hh>
#include <brick/numeric/convolve1D.hh>
#include <brick/numeric/transform2D.hh>
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
      void testLowPassEdge();
      void testVerticalEdge();

    private:

      double m_defaultTolerance;

    }; // class Iso12233Test


    /* ============== Member Function Definititions ============== */

    Iso12233Test::
    Iso12233Test()
      : brick::test::TestFixture<Iso12233Test>("Iso12233Test"),
        m_defaultTolerance(1.0E-8)
    {
      // BRICK_TEST_REGISTER_MEMBER(testLowPassEdge);
      BRICK_TEST_REGISTER_MEMBER(testVerticalEdge);
    }


    void
    Iso12233Test::
    testLowPassEdge()
    {
      constexpr std::size_t patchWidth = 128;
      constexpr std::size_t patchHeight = 100;
      constexpr std::size_t windowWidth = 64;
      constexpr double darkColor = 100.0;
      constexpr double lightColor = 200.0;
      constexpr double kernelSigma = 0.25;

      // Create a single row with a dark-to-light transition.
      Array1D<double> prototypeRow(patchWidth);
      {
        std::size_t cc = 0;
        while(cc < ((patchWidth / 2) - (windowWidth / 8))) {
          prototypeRow[cc] = darkColor;
          ++cc;
        }
        while(cc < patchWidth) {
          prototypeRow[cc] = lightColor;
          ++cc;
        }
      }

      // Low-pass filter the input row.
      Array1D<double> kernel = brick::numeric::getGaussian1D<double>(
        kernelSigma);
      kernel /= brick::numeric::sum<double>(kernel);
      Array1D<double> blurredRow = brick::numeric::convolve1D<double>(
        kernel, prototypeRow, brick::numeric::BRICK_CONVOLVE_REFLECT_SIGNAL);
      
      // Create a test image with a vertically straight up and down
      // blurred edge.  We'll rotate it later.
      Array2D<double> edgeArray(
        patchHeight, patchWidth); // Rows, columns.
      for(std::size_t rr = 0; rr < patchHeight; ++rr) {
        edgeArray.getRow(rr).copy(blurredRow);
      }

      // Create a transform that will rotate the image 15 degrees
      // around its center.
      double theta = 15.0 * brick::common::constants::radiansPerDegree;
      double cosTheta = brick::common::cosine(theta);
      double sinTheta = brick::common::sine(theta);
      double tx = ((1.0 - cosTheta) * patchWidth / 2.0
                   + sinTheta * patchHeight / 2.0);
      double ty = (-sinTheta * patchWidth / 2.0
                   + (1.0 - cosTheta) * patchHeight / 2.0);
      brick::numeric::Transform2D<double> rotation(cosTheta, -sinTheta, tx,
                                                   sinTheta, cosTheta, ty,
                                                   0.0, 0.0, 1.0);

      // Rotate the image patch.
      Image<brick::computerVision::GRAY8> rotatedImage(
        patchHeight, patchWidth); // Rows, columns.
      brick::numeric::BilinearInterpolator<double>
        interpolator(edgeArray);
      for(std::size_t rr = 0; rr < patchHeight; ++rr) {
        for(std::size_t cc = 0; cc < patchWidth; ++cc) {
          brick::numeric::Vector2D<double> outputCoord(cc, rr);
          brick::numeric::Vector2D<double> rotatedCoord =
            rotation * outputCoord;
          rotatedCoord.setValue(
            brick::common::clip(rotatedCoord.x(), 0.1, patchWidth - 1.1),
            brick::common::clip(rotatedCoord.y(), 0.1, patchHeight - 1.1));
          rotatedImage(rr, cc) = static_cast<uint8_t>(
            interpolator(rotatedCoord.y(), rotatedCoord.x()) + 0.5);
        }
      }
      
      // Try to process the image.
      Array1D<double> mtf = iso12233<double>(rotatedImage, windowWidth,
                                             [](double arg){return arg;});
    }

    
    void
    Iso12233Test::
    testVerticalEdge()
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
  brick::iso12233::Iso12233Test currentTest;
  bool result = currentTest.run();

  currentTest.testLowPassEdge();

  return (result ? 0 : 1);
}

#else

namespace {

  brick::iso12233::Iso12233Test currentTest;

}

#endif
