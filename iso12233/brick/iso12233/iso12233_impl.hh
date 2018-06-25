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

#include <brick/linearAlgebra/linearAlgebra.hh>
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
      estimateEdgeSlopeAndOffset(FloatType& slope, FloatType& offset,
                                 Array2D<FloatType> const& reflectanceImage,
                                 std::size_t windowSize);

      // Select windowSize pixels from each row of reflectanceImage,
      // with the windows centered at column coordinates given by
      // centroidArray, and return the results stacked into a
      // rectangular array.
      template <class FloatType>
      Array2D<FloatType>
      selectEdgeWindows(Array2D<FloatType> const& reflectanceImage,
                        Array1D<FloatType> const& centroidArray,
                        std::size_t windowSize);


      // Align the rows of reflectanceImage and combine them into a
      // single supersampled version of the (nearly vertical) edge.
      template <class FloatType>
      Array1D<FloatType>
      shiftAndCombineRows(Array2D<FloatType> const& reflectanceImage,
                          std::size_t const windowSize,
                          FloatType const slope,
                          FloatType const offset);

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
             ConversionFunction const& oecf)
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
      // columns).  After this call, y ~= slope*x + offset, where y is
      // column location of the edge, and x is the row number.  This
      // reverses the conventional meanings of x and y (which normally
      // mean column and row, respectively), but is consistent with
      // the variable definitions in the standard.
      FloatType slope = 0.0;
      FloatType offset = 0.0;
      privateCode::estimateEdgeSlopeAndOffset(
        slope, offset, reflectanceImage, windowSize);

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
        reflectanceImage, windowSize, slope, offset);
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
                                   / FloatType(0.5));
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

        // Finally, correct for the bias introduced by our discrete
        // three-element derivative.  The standard specifies the inverse
        // sin function used below without justification.  It looks
        // about right to me, but I haven't actually verified that this
        // is correct.
        FloatType twoPiOverN = (
          FloatType(2) * FloatType(brick::common::constants::pi)
          / FloatType(sfr.size()));
        for(std::size_t ii = 0; ii < sfr.size(); ++ii) {
          FloatType weight = 1.0 / (FloatType(ii) * twoPiOverN);
          weight = std::min(weight, FloatType(10));
          sfr[ii] *= weight;
        }

        return sfr;
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
      estimateEdgeSlopeAndOffset(FloatType& slope, FloatType& offset,
                                 Array2D<FloatType> const& reflectanceImage,
                                 std::size_t windowSize)
      {
        // To make things easier later in the function, we require that
        // reflectanceImage be wider than windowSize, and process only a
        // subset of the image data.
        std::size_t windowStartColumn =
          (reflectanceImage.columns() - windowSize) / 2;
        std::size_t windowStopColumn = windowStartColumn + windowSize;

        // Paragraph 6.2.3.2 of the standard.

        // Each line is multiplied by a Hamming window.  Note that we
        // make this bigger than windowSize because we are going to
        // convolve with an FIR derivative filter, and we want to ignore
        // the edges, where the filter extends outside windowedImage.
        if(windowSize >= (reflectanceImage.columns() - 2)) {
          BRICK_THROW(brick::common::ValueException,
                      "estimateEdgeSlopeAndOffset()",
                      "Argument reflectanceImage must have at least "
                      "(windowSize + 2) columns.");
        }
        std::size_t roiStartColumn = windowStartColumn - 1;
        std::size_t roiStopColumn = windowStopColumn + 1;
        // SubArrays do deep copies.
        Array2D<FloatType> windowedImage = brick::numeric::subArray(
          reflectanceImage, brick::numeric::Slice(),
          brick::numeric::Slice(roiStartColumn, roiStopColumn));
        Array1D<FloatType> windowFunction =
          brick::numeric::getHammingWindow1D<FloatType>(
            roiStopColumn - roiStartColumn);
        for(std::size_t rr = 0; rr < reflectanceImage.rows(); ++rr) {
          windowedImage.getRow(rr) *= windowFunction;
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

        // Estimate slope and offset of the line.  Defining rowIndices
        // to be a vector [0, 1, 2, ...]^T, we're looking for slope and
        // offset such that
        //
        //   centroidArray^T = slope * rowIndices^T + offset.
        Array1D<FloatType> rowIndices(centroidArray.size());
        for(std::size_t ii = 0; ii < rowIndices.size(); ++ii) {
          rowIndices[ii] = static_cast<FloatType>(ii);
        }
        std::pair<FloatType, FloatType> slope_offset =
          brick::linearAlgebra::linearFit(rowIndices, centroidArray);


        // Paragraph 6.2.3.3 of the standard.

        // Throw away the noise, and retain only the best fit location
        // of the line.
        centroidArray = slope_offset.first * rowIndices + slope_offset.second;

        // Repeat the Hamming window multiplication, but this time
        // center each window at the estimated line position.
        windowedImage = privateCode::selectEdgeWindows(
          reflectanceImage, centroidArray, roiStopColumn - roiStartColumn);
        for(std::size_t rr = 0; rr < reflectanceImage.rows(); ++rr) {
          windowedImage.getRow(rr) *= windowFunction;
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
        }

        // Paragraph 6.2.3.5 of the standard.
        // Re-estimate slope and offset of the line, just as above.
        slope_offset = brick::linearAlgebra::linearFit(
          rowIndices, centroidArray);

        // Pass result back to calling context.
        slope = slope_offset.first;
        offset = slope_offset.second;
      }


      // Select windowSize pixels from each row of reflectanceImage,
      // with the windows centered at column coordinates given by
      // centroidArray, and return the results stacked into a
      // rectangular array.
      template <class FloatType>
      Array2D<FloatType>
      selectEdgeWindows(Array2D<FloatType> const& reflectanceImage,
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
        }

        return outputImage;
      }


      // Align the rows of reflectanceImage and combine them into a
      // single supersampled version of the (nearly vertical) edge.
      template <class FloatType>
      Array1D<FloatType>
      shiftAndCombineRows(Array2D<FloatType> const& reflectanceImage,
                          std::size_t const windowSize,
                          FloatType const slope,
                          FloatType const offset)
      {
        // The standard requires us to supersample by a factor of four.
        // We choose to put the edge in the middle of the output row.
        std::size_t outputWindowSize = windowSize << 2;
        FloatType outputEdgeLocation = (static_cast<FloatType>(outputWindowSize)
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

          // Variable edgeLocation will be a number like "55.21"
          FloatType edgeLocation = offset + slope * rr;

          // There's a linear mapping between output row position and
          // input row position.  It has the form:
          //
          // @code
          //   outputPosition = 4 * inputPosition + delta.
          // @endcode
          //
          // The factor of 4 is required by the standard.  Now that we
          // have edge location in this row, we can solve for delta.
          FloatType delta = outputEdgeLocation - (FloatType(4.0) * edgeLocation);

          // Invert the linear equation to find out where the first and
          // last pixels of outputRow project to in inputRow.  The extra
          // 0.5 added in each row is to make the static cast round to
          // the nearest integer, rather than truncating.
          std::size_t inputStartColumn = static_cast<size_t>(
            ((FloatType(0.0) - delta) / FloatType(4.0))
            + FloatType(0.5));
          std::size_t inputStopColumn = static_cast<size_t>(
            ((FloatType(outputWindowSize - 1) - delta) / FloatType(4.0))
            + FloatType(0.5));

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

          // Adjust start and stop columns to avoid indexing issues in
          // the output row.  Remember that size_t is unsigned, so
          // negative numbers roll over to large positive numbers.
          std::size_t outputStartColumn = 0;
          while(1) {
            outputStartColumn = static_cast<std::size_t>(
              inputStartColumn * 4 + delta + 0.5);
            if(outputStartColumn < outputWindowSize) {break;}
            ++inputStartColumn;
          }
          std::size_t outputStopColumn = 0;
          while(1) {
            outputStopColumn = static_cast<std::size_t>(
              inputStopColumn * 4 + delta + 0.5);
            if(outputStopColumn < outputWindowSize) {break;}
            --inputStopColumn;
          }
          if(inputStartColumn >= inputStopColumn) {
            BRICK_THROW(brick::common::LogicException,
                        "shiftAndCombineRows()",
                        "Start and stop columns are not strictly ordered.");
          }

          // Now copy and accumulate this row's data.
          std::size_t occ = outputStartColumn;
          for(std::size_t cc = inputStartColumn; cc < inputStopColumn; ++cc) {
            outputRow[occ] += reflectanceImage(rr, cc);
            counts[occ] += 1;
            occ += 4;
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

    } // namespace privateCode

  } // namespace iso12233

} // namespace brick

#endif /* #ifndef BRICK_ISO12233_ISO12233_IMPL_HH */
