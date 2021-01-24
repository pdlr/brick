/**
***************************************************************************
* @file brick/numeric/test/fftTest.cc
*
* Source file defining FFTTest class.
*
* Copyright (C) 2017 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <complex>

#include <brick/common/functional.hh>
#include <brick/numeric/fft.hh>
#include <brick/numeric/numericTraits.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class FFTTest
      : public brick::test::TestFixture<FFTTest> {

    public:

      FFTTest();
      ~FFTTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testComputeFFT_radix2();
      void testComputeFFT_result();
      void testComputeFFT_singleFrequency();

    private:



      double m_defaultTolerance;
      double m_relaxedTolerance;

    }; // class FFTTest


    /* ============== Member Function Definititions ============== */

    FFTTest::
    FFTTest()
      : brick::test::TestFixture<FFTTest>("FFTTest"),
        m_defaultTolerance(1.0E-10),
        m_relaxedTolerance(1.0E-5)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testComputeFFT_radix2);
      BRICK_TEST_REGISTER_MEMBER(testComputeFFT_result);
      BRICK_TEST_REGISTER_MEMBER(testComputeFFT_singleFrequency);
    }


    void
    FFTTest::
    testComputeFFT_radix2()
    {
      Array1D< std::complex<double> > inputSignal(96);
      inputSignal = 0.0;
      BRICK_TEST_ASSERT_EXCEPTION(common::NotImplementedException,
                                  computeFFT(inputSignal));
    }


    void
    FFTTest::
    testComputeFFT_result()
    {
      double constexpr twoPi = brick::common::constants::twoPi;
      std::size_t constexpr signalLength = 1024;

      // Define the fourier representation of a signal.
      Array1D<double> referenceAmplitudes(signalLength);
      Array1D<double> referencePhases(signalLength);
      for(std::size_t ii = 0; ii < signalLength; ++ii) {
        referenceAmplitudes[ii] = (1.0 - double(ii) / signalLength);
        referencePhases[ii] = (double(ii) / signalLength) * twoPi;
      }

      // Compute a signal that has the fourier representation we just
      // defined.
      Array1D< std::complex<double> > inputSignal(signalLength);

      // For each element of the sequence.
      for(std::size_t ii = 0; ii < signalLength; ++ii) {
        inputSignal[ii] = std::complex<double>(0.0, 0.0);

        // For each component frequency.
        for(std::size_t jj = 0; jj < signalLength; ++jj) {
          // What's the frequency in radians per sample.
          double frequency = twoPi * double(jj) / signalLength;

          // And what are the real and imaginary parts of this component.
          double theta = ii * frequency + referencePhases[jj];
          double realPart = (referenceAmplitudes[jj] *
                             brick::common::cosine(theta));
          double imaginaryPart = (referenceAmplitudes[jj] *
                                  brick::common::sine(theta));

          // Notice that we're just Naively (O(N^2)) doing an inverse
          // DFT here, so we need a factor of 1/N to make the
          // amplitudes work out.
          realPart /= static_cast<double>(signalLength);
          imaginaryPart /= static_cast<double>(signalLength);

          // Fill in the relevant signal element.
          inputSignal[ii] += std::complex<double>(realPart, imaginaryPart);
        }
      }

      // std::cout << "\nSignal: " << inputSignal << std::endl;

      // Now do the FFT.
      Array1D< std::complex<double> > fft = computeFFT(inputSignal);

      // std::cout << "\nFFT: " << fft << std::endl;

      // Check that the result is correct.
      for(std::size_t ii = 0; ii < signalLength; ++ii) {
        double amplitude = std::abs(fft[ii]);
        double phase = brick::common::arctangent2(
          fft[ii].imag(), fft[ii].real());
        while(phase < 0.0) {phase += twoPi;}
        while(phase >= twoPi) {phase -= twoPi;}

        try {
          BRICK_TEST_ASSERT(
            approximatelyEqual(amplitude, referenceAmplitudes[ii],
                               this->m_defaultTolerance));
          if(amplitude > brick::numeric::NumericTraits<double>::epsilon()) {
            BRICK_TEST_ASSERT(
              approximatelyEqual(phase, referencePhases[ii],
                                 this->m_defaultTolerance));
          }
        } catch(...) {
          std::cout << "FFTTest::testComputeFFT_result()\n"
                    << "  Failure at element " << ii << ":\n"
                    << "  amplitude " << amplitude
                    << " at phase " << phase << " vs. reference: \n"
                    << "  amplitude " << referenceAmplitudes[ii]
                    << " at phase " << referencePhases[ii]
                    << "  differences: "
                    << amplitude - referenceAmplitudes[ii]
                    << ", " << phase - referencePhases[ii]
                    << std::endl;
          throw;
        }
      }
    }


    void
    FFTTest::
    testComputeFFT_singleFrequency()
    {
      double constexpr onePi = brick::common::constants::pi;
      double constexpr twoPi = brick::common::constants::twoPi;
      std::size_t constexpr signalLength = 1024;

      // We'll deal with real-valued signals of a single frequency,
      // where the frequency ranges from zero to Nyquist.
      for(std::size_t ii = 0; ii < signalLength / 2; /* ++ii */ ii += 2) {
        // What's the frequency in radians per sample.
        double elementAmplitude = 1.0 / signalLength;
        double frequency = twoPi * double(ii) / signalLength;
        
        // Compute a signal that has just this one frequency.  For
        // efficiency, we might move this declaration outside of the
        // enclosing for loop, but it's clearer here.
        Array1D< std::complex<double> > inputSignal(signalLength);

        // For each element of the input sequence.
        for(std::size_t jj = 0; jj < signalLength; ++jj) {
          inputSignal[jj] = std::complex<double>(
            elementAmplitude * brick::common::cosine(jj * frequency),
            0.0);
        }

        // std::cout << "\nSignal: " << inputSignal << std::endl;
      
        // Now do the FFT.
        Array1D< std::complex<double> > fft = computeFFT(inputSignal);

        // std::cout << "\nFFT: " << fft << std::endl;
      
        // Check that the result is correct.
        for(std::size_t jj = 0; jj < signalLength; ++jj) {
          double amplitude = std::abs(fft[jj]);
          double phase = brick::common::arctangent2(
            fft[jj].imag(), fft[jj].real());
          while(phase < -onePi) {phase += twoPi;}
          while(phase >= onePi) {phase -= twoPi;}

          try {
            if(jj == ii) {

              if(ii == 0 || ii == (signalLength / 2)) {
                
                // Element 0 represents DC component.  When frequency
                // is 0, all the energy goes into this element, so the
                // Fourier coefficient should be 1.0.  Similarly,
                // element signalLength / 2 represents frequency pi
                // radians per sample, which aliases to itself.  Again
                // we expect coeffient == 1.
                BRICK_TEST_ASSERT(
                  approximatelyEqual(amplitude, 1.0, this->m_defaultTolerance));
                
              } else {

                // Energy at all other frequencies gets split between
                // the Fourier coefficient at index ii and the aliased
                // frequency (above Nyquist) at (signalLength - ii).
                // These correspond to frequencies of (2 * pi * ii /
                // N) radians per sample and (2 * pi * (1 - (ii / N)))
                // radians per sample, respectively.  We expect both
                // of these coefficients to be 0.5.
                BRICK_TEST_ASSERT(
                  approximatelyEqual(amplitude, 0.5, this->m_defaultTolerance));
                
              }

              // In either case, phase should be zero, because we didn't
              // introduce any phase shift when we constructed the signal.
              BRICK_TEST_ASSERT(
                approximatelyEqual(phase, 0.0, this->m_defaultTolerance));
              
            } else if (jj == (signalLength - ii)) {

                // This is the (2 * pi * (1 - (ii / N))) radians per
                // sample component mentioned above.
                BRICK_TEST_ASSERT(
                  approximatelyEqual(amplitude, 0.5, this->m_defaultTolerance));
                BRICK_TEST_ASSERT(
                  approximatelyEqual(phase, 0.0, this->m_defaultTolerance));

            } else {

              // All other components should have zero amplitude.  Because
              // amplitude is zero (ish), phase is ill-conditioned, so we
              // ignore phase.
              BRICK_TEST_ASSERT(
                approximatelyEqual(amplitude, 0.0, this->m_defaultTolerance));
              // BRICK_TEST_ASSERT(
              //   approximatelyEqual(phase, 0.0, this->m_defaultTolerance));
              
            }
          } catch(...) {
            std::cout << "FFTTest::testComputeFFT_singleFrequency()\n"
                      << "  Failure at freqency " << ii << ", element " << jj
                      << ":\n"
                      << "  amplitude " << amplitude
                      << " at phase " << phase << std::endl;
            throw;
          }
        }
      }
    }
    
  } //  namespace numeric

} // namespace brick


#if 1

int main(int /* argc */, char** /* argv */)
{
  brick::numeric::FFTTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::FFTTest currentTest;

}

#endif
