/**
***************************************************************************
* @file brick/iso12233/iso12233.hh
*
* Header file declaring functions that implement the ISO12233 MTF
* measurement algorithms.
*
* Copyright (C) 2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_ISO12233_ISO12233_HH
#define BRICK_ISO12233_ISO12233_HH

#include <brick/iso12233/symbolImports.hh>

namespace brick {

  namespace iso12233 {


    /**
     ** This struct controls the operation of the e-SFR algorithm.
     ** The default settings give what the author considers to be good
     ** performance.  If you want things to run strictly to the
     ** standard, use Iso12233Config::getStandardCompliant() to set up
     ** an Iso12233Config instance for you.
     **/
    struct Iso12233Config {

      // Configuration variables.
      bool useInitialHammingWindow = false;
      bool useSecondHammingWindow = false;

      static Iso12233Config
      getStandardCompliant() {
        Iso12233Config result;
        result.useInitialHammingWindow = true;
        result.useSecondHammingWindow = true;
        return result;
      }
    };


    /**
     * This function implements the ISO-12233 e-SFR algorithm as
     * closely as possible.
     *
     * @param inputPatch This argument is the image patch to be
     * analyzed.  It must contain a nearly-vertical dark->light edge,
     * with the dark portion on the left side of the patch.  The patch
     * should be significantly wider (more than 2x) than argument
     * windowSize.  The patch should have plenty (say, more than 10,
     * fewer than 1000) rows.
     *
     * @param windowSize This argument must (for now) be a power of
     * two, and is the size of the window surrounding the edge that
     * should be evaluated.  Choose a number high enough to be sure
     * that any lingering effects of being near an edge have died out
     * by the time you are windowSize / 2 from the edge itself.  The
     * larger this number, the less your result will be affected by
     * FFT artifacts, so choose this to be at least 10 times the width
     * of your blurred edge.
     *
     * @param conversionFunction This function should implement the
     * Opto-Electronic Conversion Function (OECF) of the camera.  When
     * a pixel is passed to its application operator, it should return
     * a FloatType indicating the corresponding corrected intensity.
     * For example, you should be able to write "FloatType intensity =
     * oecf(inputPatch(rr, cc));" This is also a good way to select
     * specific color channels if you want to measure their MTFs
     * independently.
     *
     * @return The return value is the computed MTF of the imaging
     * system, assuming the physical edge being observed is a perfect
     * instantaneous transition from dark to light.
     */
    template <class FloatType,
              ImageFormat InputFormat,
              class ConversionFunction>
    Array1D<FloatType>
    iso12233(Image<InputFormat> const& inputPatch,
             std::size_t windowSize,
             ConversionFunction const& oecf,
             Iso12233Config const& config = Iso12233Config());

  } // namespace iso12233

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/iso12233/iso12233_impl.hh>

#endif /* #ifndef BRICK_ISO12233_ISO12233_HH */
