/**
***************************************************************************
* @file brick/iso12233/iso12233_impl.hh
*
* Header file defining functions that implement the ISO12233 MTF
* measurement algorithms.
*
* Copyright (C) 2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_ISO12233_ISO12233_IMPL_HH
#define BRICK_ISO12233_ISO12233_IMPL_HH

// This file is included by iso12233.hh, and should not be directly included
// by user code, so no need to include iso12233.hh here.
//
// #include <brick/iso12233/iso12233.hh>

#include <complex>

// #include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/computerVision/fitPolynomial.hh>
#include <brick/numeric/fft.hh>
#include <brick/numeric/numericTraits.hh>
#include <brick/numeric/sampledFunctions.hh>
#include <brick/numeric/subArray2D.hh>
#include <brick/numeric/utilities.hh>

// Locally defined templates live in the namespace "privateCode."
// Declarations go here.  Definitions are at the bottom of the file.
namespace brick {

  namespace iso12233 {

    namespace privateCode {

      // Circularly rotate the input signal so that its maximum is
      // (nearly) at the center.
      template <class FloatType>
      Array1D<FloatType>
      centerMaximum(Array1D<FloatType> const& inputRow);

      // Approximate the derivative along a row using [-1/2, 0, 1/2]
      // FIR filter.
      template <class FloatType>
      Array1D<FloatType>
      computeDerivative(Array1D<FloatType> const& inputRow);

      template <class FloatType>
      Array1D<FloatType>
      computeSFR(Array1D<FloatType> const& lineSpreadFunction);

      // Approximate the derivative along each row using [-1/2, 1/2]
      // FIR filter.
      template <class FloatType>
      Array2D<FloatType>
      computeTwoElementDerivative(Array2D<FloatType> const& inputImage);

      template <class FloatType>
      void
      estimateEdgeShape(brick::numeric::Polynomial<FloatType>& edgeShape,
                        Array2D<FloatType> const& reflectanceImage,
                        std::size_t windowSize,
                        Iso12233Config const& config);


      // Compensate for non-ideal 3-element derivative using
      // sinc-based weights.
      template <class FloatType>
      void
      reweightSFR(Array1D<FloatType>& sfr);


      // Select windowSize pixels from each row of reflectanceImage,
      // with the windows centered at column coordinates given by
      // centroidArray, and return the results stacked into a
      // rectangular array.
      template <class FloatType>
      Array2D<FloatType>
      selectEdgeWindows(Array1D<std::size_t>& pixelShiftArray,
                        Array2D<FloatType> const& reflectanceImage,
                        Array1D<FloatType> const& centroidArray,
                        std::size_t windowSize);


      // Align the rows of reflectanceImage and combine them into a
      // single supersampled version of the (nearly vertical) edge.
      template <class FloatType>
      Array1D<FloatType>
      shiftAndCombineRows(
        Array2D<FloatType> const& reflectanceImage,
        std::size_t const windowSize,
        brick::numeric::Polynomial<FloatType> const& edgeShape);

      // Subsamble the input signal by a factor of four, using a binomial
      // low-pass filter.
      template <class FloatType>
      Array1D<FloatType>
      subsampleSfr(Array1D<FloatType> const& inputSfr);

    } // namespace privateCode

  } // namespace iso12233

} // namespace brick


namespace brick {

  namespace iso12233 {

    template <class FloatType,
              ImageFormat InputFormat,
              class ConversionFunction>
    Array1D<FloatType>
    iso12233(Image<InputFormat> const& inputPatch,
             std::size_t windowSize,
             ConversionFunction const& oecf,
             Iso12233Config const& config)
    {
      // Argument checking.
      // TBD(xxx): initial check on windowSize.

      // Paragraph 6.2.2 of the standard: undo the effects of gamma
      // correction, etc.
      Array2D<FloatType> reflectanceImage(inputPatch.rows(),
                                          inputPatch.columns());
      std::transform(inputPatch.begin(), inputPatch.end(),
                     reflectanceImage.begin(), oecf);

      // Paragraph 6.2.3 of the standard.  We assume the edge is close
      // to vertical in the image (nearly aligned with the image
      // columns).  The standard calls for computing slope and offset
      // so that y ~= slope*x + offset, where y is column location of
      // the edge, and x is the row number.  This reverses the
      // conventional meanings of x and y (which normally mean column
      // and row, respectively), but is consistent with the variable
      // definitions in the standard.  By representing slope and
      // offset as a polynomial, we make it easy to use higher-order
      // shapes (quadratics, cubics) if the lens has significant
      // distortion.  The degree of this polynomial is controlled by
      // the config member variable polynomialOrder.
      brick::numeric::Polynomial<FloatType> edgeShape;
      privateCode::estimateEdgeShape(
        edgeShape, reflectanceImage, windowSize, config);

      // Paragraph 6.2.4 of the standard: align the edges in all of (or a
      // plurality of) the rows.  Because the line crosses each row at a
      // different (non-integer) column location, this alignment gives us
      // a bunch of differently-phased samplings of the edge.  The shifted
      // pixel locations of each row are then projected into new pixel
      // "bins," which are of finer pitch than the original pixel array,
      // creating a supersampled estimate of the edge shape.  The array
      // lineSpreadFunction is just the finite-differences derivative of
      // edgeSpreadFunction.
      Array1D<FloatType> edgeSpreadFunction = privateCode::shiftAndCombineRows(
        reflectanceImage, windowSize, edgeShape);
      Array1D<FloatType> lineSpreadFunction = privateCode::computeDerivative(
        edgeSpreadFunction);

      // Paragraph 6.2.5 of the standard: Compute spectral frequency response.
      Array1D<FloatType> centeredLineSpreadFunction =
        privateCode::centerMaximum(lineSpreadFunction);
      Array1D<FloatType> windowFunction =
        brick::numeric::getHammingWindow1D<FloatType>(
          centeredLineSpreadFunction.size());
      centeredLineSpreadFunction *= windowFunction;
      Array1D<FloatType> sfr = privateCode::computeSFR(
        centeredLineSpreadFunction);

      return sfr;
    }

  } // namespace iso12233

} // namespace brick


namespace brick {

  namespace iso12233 {

    namespace privateCode {

      // Circularly rotate the input signal so that its maximum is
      // (nearly) at the center.
      template <class FloatType>
      Array1D<FloatType>
      centerMaximum(Array1D<FloatType> const& inputRow)
      {
        std::size_t maxIndex = brick::numeric::argmax(inputRow);
        std::size_t newMaxIndex = inputRow.size() / 2;
        int newFromOld = (static_cast<int>(newMaxIndex)
                          - static_cast<int>(maxIndex));

        Array1D<FloatType> outputRow(inputRow.size());
        if(newFromOld > 0) {
          // Maximum is on the "left" half of inputRow.
          std::size_t newColumn = newFromOld;
          std::size_t oldColumn = 0;
          while(newColumn < outputRow.size()) {
            outputRow[newColumn] = inputRow[oldColumn];
            ++newColumn;
            ++oldColumn;
          }
          newColumn = 0;
          while(oldColumn < inputRow.size()) {
            outputRow[newColumn] = inputRow[oldColumn];
            ++newColumn;
            ++oldColumn;
          }
        } else {
          // Maximum is on the "right" half of inputRow.
          std::size_t newColumn = 0;
          std::size_t oldColumn = inputRow.size() - newFromOld;
          while(oldColumn < inputRow.size()) {
            outputRow[newColumn] = inputRow[oldColumn];
            ++newColumn;
            ++oldColumn;
          }
          oldColumn = 0;
          while(newColumn < inputRow.size()) {
            outputRow[newColumn] = inputRow[oldColumn];
            ++newColumn;
            ++oldColumn;
          }
        }
        return outputRow;
      }


      // Approximate the derivative along a row using [-1/2, 0, 1/2]
      // FIR filter.
      template <class FloatType>
      Array1D<FloatType>
      computeDerivative(Array1D<FloatType> const& inputRow)
      {
        // Convolve with the derivative kernel.
        Array1D<FloatType> derivativeRow(inputRow.size());
        for(std::size_t column = 1; column < inputRow.size() - 1; ++column) {
          derivativeRow[column] = ((inputRow[column + 1] - inputRow[column - 1])
                                   * FloatType(0.5));
        }

        // Copy first and last elements to fill out the array.
        derivativeRow[0] = derivativeRow[1];
        derivativeRow[inputRow.size() - 1] = derivativeRow[inputRow.size() - 2];
        return derivativeRow;
      }


      template <class FloatType>
      Array1D<FloatType>
      computeSFR(Array1D<FloatType> const& lineSpreadFunction)
      {
        // Our FFT routine only knows about complex numbers.
        Array1D<std::complex<FloatType> > complexLSF(lineSpreadFunction.size());
        std::transform(lineSpreadFunction.begin(), lineSpreadFunction.end(),
                       complexLSF.begin(),
                       [](FloatType const& arg0) {
                         return std::complex<FloatType>(arg0, FloatType(0));
                       });

        // Here we compute the Discrete Fourier Transform (DFT).
        Array1D<std::complex<FloatType> > dft =
          brick::numeric::computeFFT(complexLSF);

        // Spacial Frequency Response (SFR) is the normalized modulus of
        // the DFT.
        Array1D<FloatType> sfr(dft.size());
        std::transform(dft.begin(), dft.end(), sfr.begin(),
                       [](std::complex<FloatType> const& arg0) {
                         return std::abs(arg0);
                       });

        // Now normalize so that DC component is 1.0.
        FloatType epsilon = brick::numeric::NumericTraits<FloatType>::epsilon();
        if(brick::numeric::absoluteValue(sfr[0]) <= epsilon) {
          BRICK_THROW(brick::common::ValueException,
                      "computeSFR()",
                      "Spacial frequency response appears to have no "
                      "DC component.");
        }
        sfr /= sfr[0];

        // correct for the bias introduced by our discrete three-element
        // derivative.
        reweightSFR(sfr);

        // Finally, return the low-frequency component.
        return subsampleSfr(sfr);
      }


      // Approximate the derivative along each row using [-1/2, 1/2]
      // FIR filter.
      template <class FloatType>
      Array2D<FloatType>
      computeTwoElementDerivative(Array2D<FloatType> const& inputImage)
      {
        Array2D<FloatType> derivativeImage(inputImage.rows(),
                                           inputImage.columns());
        for(std::size_t row = 0; row < inputImage.rows(); ++row) {
          derivativeImage(row, 0) = FloatType(0.0);
          for(std::size_t column = 1; column < inputImage.columns();
              ++column) {
            derivativeImage(row, column) =
              ((inputImage(row, column) - inputImage(row, (column - 1)))
               / FloatType(2.0));
          }
          derivativeImage(row, 0) = derivativeImage(row, 1);
        }
        return derivativeImage;
      }


      template <class FloatType>
      void
      estimateEdgeShape(brick::numeric::Polynomial<FloatType>& edgeShape,
                        Array2D<FloatType> const& reflectanceImage,
                        std::size_t windowSize,
                        Iso12233Config const& config)
      {
        // To make things easier later in the function, we require that
        // reflectanceImage be wider than windowSize, and process only a
        // subset of the image data.
        if(windowSize >= (reflectanceImage.columns() - 2)) {
          BRICK_THROW(brick::common::ValueException,
                      "estimateEdgeShape()",
                      "Argument reflectanceImage must have at least "
                      "(windowSize + 2) columns.");
        }
        // std::size_t windowStartColumn =
        //   (reflectanceImage.columns() - windowSize) / 2;
        // std::size_t windowStopColumn = windowStartColumn + windowSize;

        // // xxx
        // // Note that we
        // // make this bigger than windowSize because we are going to
        // // convolve with an FIR derivative filter, and we want to ignore
        // // the edges, where the filter extends outside windowedImage.
        // // SubArrays do deep copies.
        // std::size_t roiStartColumn = windowStartColumn - 1;
        // std::size_t roiStopColumn = windowStopColumn + 1;
        // Array2D<FloatType> windowedImage = brick::numeric::subArray(
        //   reflectanceImage, brick::numeric::Slice(),
        //   brick::numeric::Slice(roiStartColumn, roiStopColumn));

        // Paragraph 6.2.3.2 of the standard.
        // Each line is multiplied by a Hamming window.
        Array2D<FloatType> windowedImage = reflectanceImage.copy();
        if(config.useInitialHammingWindow) {
          Array1D<FloatType> windowFunction =
            brick::numeric::getHammingWindow1D<FloatType>(
              windowedImage.columns());
          for(std::size_t rr = 0; rr < reflectanceImage.rows(); ++rr) {
            windowedImage.getRow(rr) *= windowFunction;
          }
        }

        // Approximate the derivative along each row using [-1/2, 1/2]
        // FIR filter.
        Array2D<FloatType> derivativeImage = computeTwoElementDerivative(
          windowedImage);

        // The 1D centroid is computed for each line.
        Array1D<FloatType> centroidArray(derivativeImage.rows());
        for(std::size_t rr = 0; rr < derivativeImage.rows(); ++rr) {
          centroidArray[rr] = brick::numeric::getCentroid<FloatType>(
            derivativeImage.getRow(rr));

          // Compensate for the fact that our derivative filter had an
          // even number of elements.  This is not described in the
          // body of the standard, but ought to be. It _is_ reflected in
          // appendix D of the standard, though.
          centroidArray[rr] -= FloatType(0.5);
        }

        // The standard calls for us to estimate slope and offset of
        // the line here.  Defining rowIndices to be a vector
        // [0, 1, 2, ...]^T, we would be looking for slope and offset
        // such that
        //
        //   centroidArray^T = slope * rowIndices^T + offset.
        //
        // In practice, we're normally working with cameras that have
        // some degree of lens distortion, so we allow the user to
        // specify a polynomial of arbitrary order.  Setting
        // polynomialOrder to 1 gives the (linear) fit called for by
        // the standard and described above in this comment.  Setting
        // polynomialOrder to 2 gives a quadratic fit, which is likely
        // flexible enough for most applications.
        Array1D<FloatType> rowIndices(centroidArray.size());
        for(std::size_t ii = 0; ii < rowIndices.size(); ++ii) {
          rowIndices[ii] = static_cast<FloatType>(ii);
        }
        // std::pair<FloatType, FloatType> slope_offset =
        //   brick::linearAlgebra::linearFit(rowIndices, centroidArray);
        edgeShape = brick::computerVision::fitPolynomial<FloatType>(
          rowIndices.begin(), rowIndices.end(), centroidArray.begin(),
          config.polynomialOrder);

        // Paragraph 6.2.3.3 of the standard.

        // Throw away the noise, and retain only the best fit location
        // of the line.
        //
        // centroidArray = (slope_offset.first * rowIndices
        //                  + slope_offset.second);
        //
        // This call just applies the polynomial to each row index,
        // and stores the result in centroidArray.
        std::transform(rowIndices.begin(), rowIndices.end(),
                       centroidArray.begin(), edgeShape);

        // Repeat the Hamming window multiplication, but this time
        // center each window at the estimated line position, and
        // restrict attention to the user specified window size.  Note
        // that we grab a little more image data than just windowSize
        // because we are going to convolve with an FIR derivative
        // filter, and we want to ignore the edges, where the filter
        // extends outside windowSize.  Argument pixelShiftArray tells
        // us which pixels got selected for each row of windowedImage,
        // so we can correctly calculate the line polynomial later.
        Array1D<std::size_t> pixelShiftArray;
        windowedImage = privateCode::selectEdgeWindows(
          pixelShiftArray, reflectanceImage, centroidArray, windowSize + 2);
        if(config.useSecondHammingWindow) {
          Array1D<FloatType> windowFunction =
            brick::numeric::getHammingWindow1D<FloatType>(
              windowedImage.columns());
          for(std::size_t rr = 0; rr < reflectanceImage.rows(); ++rr) {
            windowedImage.getRow(rr) *= windowFunction;
          }
        }

        // Approximate the derivative along each row again using
        // [-1/2, 1/2] FIR filter.
        derivativeImage = computeTwoElementDerivative(windowedImage);

        // Paragraph 6.2.3.4 of the standard.
        // The 1D centroid is computed for each line.
        for(std::size_t rr = 0; rr < derivativeImage.rows(); ++rr) {
          centroidArray[rr] = brick::numeric::getCentroid<FloatType>(
            derivativeImage.getRow(rr));

          // Compensate for the fact that our derivative filter had an
          // even number of elements.  This is not described in the
          // body of the standard, but ought to be. It _is_ reflected in
          // appendix D of the standard, though.
          centroidArray[rr] -= FloatType(0.5);

          // Compensate for the fact that each row of this image was
          // shifted to center the Hamming window.
          centroidArray[rr] += static_cast<FloatType>(pixelShiftArray[rr]);
        }

        // Paragraph 6.2.3.5 of the standard.
        // Re-estimate slope and offset of the line, just as above.
        //
        // slope_offset = brick::linearAlgebra::linearFit(
        // rowIndices, centroidArray);
        //
        // As before, allow polynomial fits of arbitrary order.
        edgeShape = brick::computerVision::fitPolynomial<FloatType>(
          rowIndices.begin(), rowIndices.end(), centroidArray.begin(),
          config.polynomialOrder);
      }


      // Compensate for non-ideal 3-element derivative using
      // sinc-based weights.
      template <class FloatType>
      void
      reweightSFR(Array1D<FloatType>& sfr)
      {
        // If we were working in continuous signals, rather than
        // discrete, then taking the ideal derivative of the edge
        // image to get the line spread function would look like
        // multiplying the DFT by a linear function with zero at the
        // origin and slope == 1.0.
        //
        // Back in the real world, we convolved with [-0.5, 0, 0.5].
        // If we think of this as a sequence, h[n], of N values that is
        // nonzero only in the first and third elements, the DFT of
        // this signal becomes clear.
        //
        // The response of the derivative filter is:
        // @code
        //   H[k]   = sum_n(h[n] * exp[-j*2*pi*k*n/N])
        //          = 0.5*exp[-j*4*pi*k/N] - 0.5*exp[0]
        //          = (0.5*cos(4*pi*k/N) - 0.5) - j*(0.5*sin(4*pi*k/N))
        //   |H[k]| = sqrt(H[k] * conjugate(H[k]))
        //          = sqrt(0.25*cos(4*pi*k/N)^2 - 0.5*cos(4*pi*k/N) + 0.25)
        //                 + 0.25*sin(4*pi*k/N)^2)
        //          = 0.5 * sqrt(cos(4*pi*k/N)^2 - 2*cos(4*pi*k/N) + 1)
        //                       + sin(4*pi*k/N)^2)
        //          = 0.5 * sqrt((cos(4*pi*k/N)^2 + sin(4*pi*k/N)^2) + 1
        //                        - 2*cos(4*pi*k/N))
        //          = 0.5 * sqrt(2 - 2*cos(4*pi*k/N))
        // @endcode
        // But, we have the trig identity sin^2(theta) = (1 - cos(2*theta))/2.
        // Substituting this into the above, we have:
        // @code
        //   |H[k]| = 0.5 * sqrt(4 * (1 - cos(4*pi*k/N)) / 2)
        //          = sqrt((1 - cos(4*pi*k/N)) / 2)
        //          = sqrt(sin(2*pi*k/N)^2)
        //          = |sin(2*pi*k/N)|
        // @endcode
        //
        // To get back to the ideal (linear) differentiation, we need
        // to divide by this response, and multiply by the desired
        // response.  As discussed above, the desired response is
        // |F[k]| = 2*pi*k/N, but antisymmetric around the origin
        // because our FFT implementation puts negative frequencies in
        // the second half of the output vector.  At the origin, we
        // use L'Hopitals rule to avoid dividing by zero, and set the
        // corrective weight to 1.0.  That is, we leave sfr[0]
        // untouched.
        //
        // At sfr[sfr.size() / 2], the actual response of the
        // three-element derivative filter goes to zero, so we leave this
        // element untouched, too.
        FloatType twoPiOverN = (
          FloatType(2) * FloatType(brick::common::constants::pi)
          / FloatType(sfr.size()));
        for(std::size_t kk = 1; kk < sfr.size() / 2; ++kk) {
          FloatType kkf(kk);
          FloatType desiredResponse = kkf * twoPiOverN;
          FloatType actualResponse =
            brick::numeric::absoluteValue(
              brick::numeric::sine(kkf * twoPiOverN));
          FloatType weight = desiredResponse / actualResponse;
          sfr[kk] *= weight;
          sfr[sfr.size() - kk] *= weight;
        }
      }


      // Select windowSize pixels from each row of reflectanceImage,
      // with the windows centered at column coordinates given by
      // centroidArray, and return the results stacked into a
      // rectangular array.
      template <class FloatType>
      Array2D<FloatType>
      selectEdgeWindows(Array1D<std::size_t>& pixelShiftArray,
                        Array2D<FloatType> const& reflectanceImage,
                        Array1D<FloatType> const& centroidArray,
                        std::size_t windowSize)
      {
        std::size_t const windowSizeOverTwo = windowSize / 2;

        // Check arguments.
        if(centroidArray.size() != reflectanceImage.rows()) {
          BRICK_THROW(brick::common::ValueException,
                      "selectEdgeWindows()",
                      "Argument centroidArray must have exactly one entry "
                      "for each row of argument reflectanceImage");
        }

        // Alocate space for function output.
        Array2D<FloatType> outputImage(centroidArray.size(), windowSize);
        pixelShiftArray.reinit(centroidArray.size());

        // Iterate over each row of reflectanceImage.
        for(std::size_t rr = 0; rr < centroidArray.size(); ++rr) {

          // Figure out which pixels we're interested in.
          std::size_t const centerColumn = static_cast<std::size_t>(
            centroidArray[rr] + 0.5);
          std::size_t const startColumn = centerColumn - windowSizeOverTwo;
          std::size_t const stopColumn = startColumn + windowSize;

          // Make sure these pixels are in-bounds.
          if(centerColumn < windowSizeOverTwo
             || stopColumn >= reflectanceImage.columns()) {
            BRICK_THROW(brick::common::ValueException,
                        "selectEdgeWindows()",
                        "Edge is too close to side of input patch.  "
                        "Consider using a smaller window size or a larger "
                        "input image patch.");
          }

          // Copy the relevant pixels from this row into the output image.
          std::copy(reflectanceImage.data(rr, startColumn),
                    reflectanceImage.data(rr, stopColumn),
                    outputImage.data(rr, 0));
          pixelShiftArray[rr] = startColumn;
        }

        return outputImage;
      }


      // Align the rows of reflectanceImage and combine them into a
      // single supersampled version of the (nearly vertical) edge.
      template <class FloatType>
      Array1D<FloatType>
      shiftAndCombineRows(
        Array2D<FloatType> const& reflectanceImage,
        std::size_t const windowSize,
        brick::numeric::Polynomial<FloatType> const& edgeShape)
      {
        // The standard requires us to supersample by a factor of four.
        // We choose to put the edge in the middle of the output row.
        std::size_t outputWindowSize = windowSize << 2;
        FloatType outputEdgeLocation =
          (static_cast<FloatType>(outputWindowSize - 1)
           / FloatType(2.0));

        // We'll accumulate pixel values in these arrays, which we
        // initialize to zero.
        Array1D<FloatType> outputRow(outputWindowSize);
        Array1D<std::size_t> counts(outputWindowSize);
        outputRow = FloatType(0.0);
        counts = 0;

        // Now iterate over each input row, copying the pixel values
        // into the accumulater arrays we just created.
        for(std::size_t rr = 0; rr < reflectanceImage.rows(); ++rr) {

          // Variable edgeLocation will be a number in input pixel column
          // coordinates.  Something like "55.21"
          FloatType edgeLocation = edgeShape(static_cast<FloatType>(rr));

          // There's a linear mapping between output column coordinate and
          // input column coordinate.  It has the form:
          //
          // @code
          //   outputPosition = 4 * inputPosition + delta.
          // @endcode
          //
          // The factor of 4 is required by the standard.  Now that we
          // have edge location in this row, we can solve for delta.
          FloatType delta = (outputEdgeLocation
                             - (FloatType(4.0) * edgeLocation));

          // Invert the linear equation to find out where the first and
          // last pixels of outputRow project to in inputRow.
          std::size_t inputStartColumn = static_cast<size_t>(
            std::round((FloatType(0.0) - delta) / FloatType(4.0)));
          std::size_t inputStopColumn = static_cast<size_t>(
            std::round((FloatType(outputWindowSize - 1) - delta)
                       / FloatType(4.0)));

          // Make sure no indexing issues in the input patch.  Remember
          // that size_t is unsigned, so negative numbers roll over to
          // positive numbers.
          if(inputStopColumn >= reflectanceImage.columns()
             || inputStartColumn >= reflectanceImage.columns()) {

            BRICK_THROW(brick::common::ValueException,
                        "shiftAndCombineRows()",
                        "Edge is too close to side of input patch.  "
                        "Consider using a smaller window size or a larger "
                        "input image patch.");
          }

          // We've found the input pixel to which the first output
          // pixel most nearly maps, but remember that there are 4x as
          // many output pixels as input pixels.  The center of this
          // input pixel does not necessarily map to the first output
          // pixel.  We project in the other direction now to find out
          // which output pixel gets the contents of this first pixel.
          // If the result is out of bounds, we increment the input
          // pixel index.  Remember that size_t is unsigned, so
          // negative numbers roll over to large positive numbers.
          std::size_t outputStartColumn = 0;
          while(1) {
            outputStartColumn = static_cast<std::size_t>(
              std::round(inputStartColumn * 4 + delta));
            if(outputStartColumn < outputWindowSize) {break;}
            ++inputStartColumn;
          }

          // Now copy and accumulate this row's data.
          std::size_t occ = outputStartColumn;
          std::size_t cc = inputStartColumn;
          while(cc < reflectanceImage.columns()
                && occ < outputRow.size()) {
            outputRow[occ] += reflectanceImage(rr, cc);
            counts[occ] += 1;
            occ += 4;
            ++cc;
          }
        }

        // Take mean value of each bin in outputRow.
        for(std::size_t occ = 0; occ < outputWindowSize; ++occ) {
          if(counts[occ] == 0) {
            BRICK_THROW(brick::common::ValueException,
                        "shiftAndCombineRows()",
                        "Output row is sparse.  Please try a different edge "
                        "angle.");
          }
          outputRow[occ] /= counts[occ];
        }

        return outputRow;
      }


      // Subsamble the input signal by a factor of four, using a binomial
      // low-pass filter.
      template <class FloatType>
      Array1D<FloatType>
      subsampleSfr(Array1D<FloatType> const& inputSfr)
      {
        Array1D<FloatType> outputSfr(inputSfr.size() / 4);
        std::copy(inputSfr.begin(), inputSfr.begin() + outputSfr.size(),
                  outputSfr.begin());
        outputSfr /= outputSfr[0];
        return outputSfr;
      }

    } // namespace privateCode

  } // namespace iso12233

} // namespace brick

#endif /* #ifndef BRICK_ISO12233_ISO12233_IMPL_HH */
