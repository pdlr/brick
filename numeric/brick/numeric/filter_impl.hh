/**
***************************************************************************
* @file brick/numeric/filter_impl.hh
*
* Header file defining functions for performing filtering of arrays.
*
* Copyright (C) 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_FILTER_IMPL_HH
#define BRICK_NUMERIC_FILTER_IMPL_HH

// This file is included by filter.hh, and should not be directly included
// by user code, so no need to include filter.hh here.
//
// #include <brick/numeric/filter.hh>


namespace brick {

  namespace numeric {

    namespace privateCode {

      template<class InputType, class OutputType, class KernelType>
      void
      filterColumns_3(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_3()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 3) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_3()",
                      "Argument kernel must have 3 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startRow = kernel.size() / 2;
        unsigned int const stopRow = inputArray.rows() - startRow;
        unsigned int const numColumns = inputArray.columns();
        unsigned int const rowStep      = inputArray.getRowStep();

        // Convolve, but do it row-wise to minimize cache misses.
        for(unsigned int row = startRow; row < stopRow; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = &(outputArray(row, 0));
          for(unsigned int column = 0; column < numColumns;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (
              static_cast<KernelType>(*(inPtr - rowStep))   * kernel[0]
              + static_cast<KernelType>(*(inPtr))           * kernel[1]
              + static_cast<KernelType>(*(inPtr + rowStep)) * kernel[2]);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterColumns_5(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_5()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 5) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_5()",
                      "Argument kernel must have 5 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startRow = kernel.size() / 2;
        unsigned int const stopRow = inputArray.rows() - startRow;
        unsigned int const numColumns = inputArray.columns();
        unsigned int const rowStep      = inputArray.getRowStep();
        unsigned int const twoRowStep   = 2 * rowStep;

        // Convolve, but do it row-wise to minimize cache misses.
        for(unsigned int row = startRow; row < stopRow; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = &(outputArray(row, 0));
          for(unsigned int column = 0; column < numColumns;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (
              static_cast<KernelType>(*(inPtr - twoRowStep))   * kernel[0]
              + static_cast<KernelType>(*(inPtr - rowStep))    * kernel[1]
              + static_cast<KernelType>(*(inPtr))              * kernel[2]
              + static_cast<KernelType>(*(inPtr + rowStep))    * kernel[3]
              + static_cast<KernelType>(*(inPtr + twoRowStep)) * kernel[4]);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterColumns_7(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_7()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 7) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_7()",
                      "Argument kernel must have 7 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startRow = kernel.size() / 2;
        unsigned int const stopRow = inputArray.rows() - startRow;
        unsigned int const numColumns = inputArray.columns();
        unsigned int const rowStep      = inputArray.getRowStep();
        unsigned int const twoRowStep   = 2 * rowStep;
        unsigned int const threeRowStep = 3 * rowStep;

        // Convolve, but do it row-wise to minimize cache misses.
        for(unsigned int row = startRow; row < stopRow; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = &(outputArray(row, 0));
          for(unsigned int column = 0; column < numColumns;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (
              static_cast<KernelType>(*(inPtr - threeRowStep))   * kernel[0]
              + static_cast<KernelType>(*(inPtr - twoRowStep))   * kernel[1]
              + static_cast<KernelType>(*(inPtr - rowStep))      * kernel[2]
              + static_cast<KernelType>(*(inPtr))                * kernel[3]
              + static_cast<KernelType>(*(inPtr + rowStep))      * kernel[4]
              + static_cast<KernelType>(*(inPtr + twoRowStep))   * kernel[5]
              + static_cast<KernelType>(*(inPtr + threeRowStep)) * kernel[6]);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterColumns_9(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_9()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 9) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_9()",
                      "Argument kernel must have 9 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startRow = kernel.size() / 2;
        unsigned int const stopRow = inputArray.rows() - startRow;
        unsigned int const numColumns = inputArray.columns();
        unsigned int const rowStep      = inputArray.getRowStep();
        unsigned int const twoRowStep   = 2 * rowStep;
        unsigned int const threeRowStep = 3 * rowStep;
        unsigned int const fourRowStep  = 4 * rowStep;

        // Convolve, but do it row-wise to minimize cache misses.
        for(unsigned int row = startRow; row < stopRow; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = &(outputArray(row, 0));
          for(unsigned int column = 0; column < numColumns;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (
              static_cast<KernelType>(*(inPtr - fourRowStep))    * kernel[0]
              + static_cast<KernelType>(*(inPtr - threeRowStep)) * kernel[1]
              + static_cast<KernelType>(*(inPtr - twoRowStep))   * kernel[2]
              + static_cast<KernelType>(*(inPtr - rowStep))      * kernel[3]
              + static_cast<KernelType>(*(inPtr))                * kernel[4]
              + static_cast<KernelType>(*(inPtr + rowStep))      * kernel[5]
              + static_cast<KernelType>(*(inPtr + twoRowStep))   * kernel[6]
              + static_cast<KernelType>(*(inPtr + threeRowStep)) * kernel[7]
              + static_cast<KernelType>(*(inPtr + fourRowStep))  * kernel[8]);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterColumns_11(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_11()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 11) {
          BRICK_THROW(brick::common::ValueException, "filterColumns_11()",
                      "Argument kernel must have 11 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startRow = kernel.size() / 2;
        unsigned int const stopRow = inputArray.rows() - startRow;
        unsigned int const numColumns = inputArray.columns();
        unsigned int const rowStep      = inputArray.getRowStep();
        unsigned int const twoRowStep   = 2 * rowStep;
        unsigned int const threeRowStep = 3 * rowStep;
        unsigned int const fourRowStep  = 4 * rowStep;
        unsigned int const fiveRowStep  = 5 * rowStep;

        // Convolve, but do it row-wise to minimize cache misses.
        for(unsigned int row = startRow; row < stopRow; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = &(outputArray(row, 0));
          for(unsigned int column = 0; column < numColumns;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (
              static_cast<KernelType>(*(inPtr - fiveRowStep))  * kernel[0]
              + static_cast<KernelType>(*(inPtr - fourRowStep))  * kernel[1]
              + static_cast<KernelType>(*(inPtr - threeRowStep)) * kernel[2]
              + static_cast<KernelType>(*(inPtr - twoRowStep))   * kernel[3]
              + static_cast<KernelType>(*(inPtr - rowStep))      * kernel[4]
              + static_cast<KernelType>(*(inPtr))                * kernel[5]
              + static_cast<KernelType>(*(inPtr + rowStep))      * kernel[6]
              + static_cast<KernelType>(*(inPtr + twoRowStep))   * kernel[7]
              + static_cast<KernelType>(*(inPtr + threeRowStep)) * kernel[8]
              + static_cast<KernelType>(*(inPtr + fourRowStep))  * kernel[9]
              + static_cast<KernelType>(*(inPtr + fiveRowStep))  * kernel[10]);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterColumns_generic(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterColumns()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() % 2 != 1) {
          BRICK_THROW(brick::common::ValueException, "filterColumns()",
                      "Argument kernel must have an odd number of elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startRow = kernel.size() / 2;
        unsigned int const stopRow = inputArray.rows() - startRow;
        unsigned int const numColumns = inputArray.columns();
        unsigned int const inRowStep = inputArray.getRowStep();
        unsigned int const kernelSize = kernel.size();

        // Convolve, but do it row-wise to minimize cache misses.
        for(unsigned int row = startRow; row < stopRow; ++row) {
          InputType const* inPtr = inputArray.data(row - startRow, 0);
          OutputType* outPtr = &(outputArray(row, 0));
          for(unsigned int column = 0; column < numColumns;
              ++column, ++inPtr, ++outPtr) {
            KernelType result = KernelType(0);
            InputType const* workingPtr = inPtr;
            for(unsigned int ii = 0; ii < kernelSize; ++ii) {
              result += kernel[ii] * static_cast<KernelType>(*workingPtr);
              workingPtr += inRowStep;
            }
            *outPtr = static_cast<OutputType>(result);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterRows_3(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterRows_3()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 3) {
          BRICK_THROW(brick::common::ValueException, "filterRows_3()",
                      "Argument kernel must have 3 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startColumn = kernel.size() / 2;
        unsigned int const stopColumn = inputArray.columns() - startColumn;
        unsigned int const numRows = inputArray.rows();

        // Convolve.
        for(unsigned int row = 0; row < numRows; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = outputArray.data(row, startColumn);
          for(unsigned int column = startColumn; column < stopColumn;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (static_cast<KernelType>(inPtr[0]) * kernel[0]
                       + static_cast<KernelType>(inPtr[1]) * kernel[1]
                       + static_cast<KernelType>(inPtr[2]) * kernel[2]);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterRows_5(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterRows_5()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 5) {
          BRICK_THROW(brick::common::ValueException, "filterRows_5()",
                      "Argument kernel must have 5 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startColumn = kernel.size() / 2;
        unsigned int const stopColumn = inputArray.columns() - startColumn;
        unsigned int const numRows = inputArray.rows();

        // Convolve.
        for(unsigned int row = 0; row < numRows; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = outputArray.data(row, startColumn);
          for(unsigned int column = startColumn; column < stopColumn;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (static_cast<KernelType>(inPtr[0]) * kernel[0]
                       + static_cast<KernelType>(inPtr[1]) * kernel[1]
                       + static_cast<KernelType>(inPtr[2]) * kernel[2]
                       + static_cast<KernelType>(inPtr[3]) * kernel[3]
                       + static_cast<KernelType>(inPtr[4]) * kernel[4]);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterRows_7(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterRows_7()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 7) {
          BRICK_THROW(brick::common::ValueException, "filterRows_7()",
                      "Argument kernel must have 7 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startColumn = kernel.size() / 2;
        unsigned int const stopColumn = inputArray.columns() - startColumn;
        unsigned int const numRows = inputArray.rows();

        // Convolve.
        for(unsigned int row = 0; row < numRows; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = outputArray.data(row, startColumn);
          for(unsigned int column = startColumn; column < stopColumn;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (static_cast<KernelType>(inPtr[0]) * kernel[0]
                       + static_cast<KernelType>(inPtr[1]) * kernel[1]
                       + static_cast<KernelType>(inPtr[2]) * kernel[2]
                       + static_cast<KernelType>(inPtr[3]) * kernel[3]
                       + static_cast<KernelType>(inPtr[4]) * kernel[4]
                       + static_cast<KernelType>(inPtr[5]) * kernel[5]
                       + static_cast<KernelType>(inPtr[6]) * kernel[6]);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterRows_9(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterRows_9()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 9) {
          BRICK_THROW(brick::common::ValueException, "filterRows_9()",
                      "Argument kernel must have 9 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startColumn = kernel.size() / 2;
        unsigned int const stopColumn = inputArray.columns() - startColumn;
        unsigned int const numRows = inputArray.rows();

        // Convolve.
        for(unsigned int row = 0; row < numRows; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = outputArray.data(row, startColumn);
          for(unsigned int column = startColumn; column < stopColumn;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (static_cast<KernelType>(inPtr[0]) * kernel[0]
                       + static_cast<KernelType>(inPtr[1]) * kernel[1]
                       + static_cast<KernelType>(inPtr[2]) * kernel[2]
                       + static_cast<KernelType>(inPtr[3]) * kernel[3]
                       + static_cast<KernelType>(inPtr[4]) * kernel[4]
                       + static_cast<KernelType>(inPtr[5]) * kernel[5]
                       + static_cast<KernelType>(inPtr[6]) * kernel[6]
                       + static_cast<KernelType>(inPtr[7]) * kernel[7]
                       + static_cast<KernelType>(inPtr[8]) * kernel[8]);
          }
        }
      }


      template<class InputType, class OutputType, class KernelType>
      void
      filterRows_11(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Make sure outputArray is appropriately sized.
        if(outputArray.rows() != inputArray.rows()
           || outputArray.columns() != inputArray.columns()) {
          BRICK_THROW(brick::common::ValueException, "filterRows_11()",
                      "Arguments outputArray and inputArray must have "
                      "the same shape.");
        }

        // Make sure kernel is appropriately sized.
        if(kernel.size() != 11) {
          BRICK_THROW(brick::common::ValueException, "filterRows_11()",
                      "Argument kernel must have 11 elements.");
        }

        // Get all our dimensions straight.
        unsigned int const startColumn = kernel.size() / 2;
        unsigned int const stopColumn = inputArray.columns() - startColumn;
        unsigned int const numRows = inputArray.rows();

        // Convolve.
        for(unsigned int row = 0; row < numRows; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = outputArray.data(row, startColumn);
          for(unsigned int column = startColumn; column < stopColumn;
              ++column, ++inPtr, ++outPtr) {
            *outPtr = (static_cast<KernelType>(inPtr[0]) * kernel[0]
                       + static_cast<KernelType>(inPtr[1]) * kernel[1]
                       + static_cast<KernelType>(inPtr[2]) * kernel[2]
                       + static_cast<KernelType>(inPtr[3]) * kernel[3]
                       + static_cast<KernelType>(inPtr[4]) * kernel[4]
                       + static_cast<KernelType>(inPtr[5]) * kernel[5]
                       + static_cast<KernelType>(inPtr[6]) * kernel[6]
                       + static_cast<KernelType>(inPtr[7]) * kernel[7]
                       + static_cast<KernelType>(inPtr[8]) * kernel[8]
                       + static_cast<KernelType>(inPtr[9]) * kernel[9]
                       + static_cast<KernelType>(inPtr[10]) * kernel[10]);
          }
        }
      }

      template<class InputType, class OutputType, class KernelType>
      void
      filterRows_generic(
        Array2D<OutputType>& outputArray,
        Array2D<InputType> const& inputArray,
        Array1D<KernelType> const& kernel)
      {
        // Get all our dimensions straight.
        unsigned int const startColumn = kernel.size() / 2;
        unsigned int const stopColumn = inputArray.columns() - startColumn;
        unsigned int const numRows = inputArray.rows();
        unsigned int const kernelSize = kernel.size();

        // Convolve.
        for(unsigned int row = 0; row < numRows; ++row) {
          InputType const* inPtr = inputArray.data(row, 0);
          OutputType* outPtr = outputArray.data(row, startColumn);
          for(unsigned int column = startColumn; column < stopColumn;
              ++column, ++inPtr, ++outPtr) {
            KernelType result = KernelType(0);
            for(unsigned int ii = 0; ii < kernelSize; ++ii) {
              result += kernel[ii] * static_cast<KernelType>(inPtr[ii]);
            }
            *outPtr = static_cast<OutputType>(result);
          }
        }
      }

    } // namespace privateCode

    template<class InputType, class OutputType, class KernelType>
    void
    filterColumns(
      Array2D<OutputType>& outputArray,
      Array2D<InputType> const& inputArray,
      Array1D<KernelType> const& kernel)
    {
      // Make sure outputArray is appropriately sized.
      if(outputArray.rows() != inputArray.rows()
         || outputArray.columns() != inputArray.columns()) {
        BRICK_THROW(brick::common::ValueException, "filterColumns()",
                    "Arguments outputArray and inputArray must have "
                    "the same shape.");
      }

      // Make sure kernel is appropriately sized.
      if(kernel.size() % 2 != 1) {
        BRICK_THROW(brick::common::ValueException, "filterColumns()",
                    "Argument kernel must have an odd number of elements.");
      }

      switch(kernel.size()) {
      case 3:
        privateCode::filterColumns_3(outputArray, inputArray, kernel);
        break;
      case 5:
        privateCode::filterColumns_5(outputArray, inputArray, kernel);
        break;
      case 7:
        // Fall through to faster generic routine.
        //
        // privateCode::filterColumns_7(outputArray, inputArray, kernel);
        // break;
      case 9:
        // Fall through to faster generic routine.
        //
        // privateCode::filterColumns_9(outputArray, inputArray, kernel);
        // break;
      case 11:
        // Fall through to faster generic routine.
        //
        // privateCode::filterColumns_11(outputArray, inputArray, kernel);
        // break;
      default:
        privateCode::filterColumns_generic(outputArray, inputArray, kernel);
        break;
      }
    }


    template<class InputType, class OutputType, class KernelType>
    void
    filterRows(
      Array2D<OutputType>& outputArray,
      Array2D<InputType> const& inputArray,
      Array1D<KernelType> const& kernel)
    {
      // Make sure outputArray is appropriately sized.
      if(outputArray.rows() != inputArray.rows()
         || outputArray.columns() != inputArray.columns()) {
        BRICK_THROW(brick::common::ValueException, "filterRows()",
                    "Arguments outputArray and inputArray must have "
                    "the same shape.");
      }

      // Make sure kernel is appropriately sized.
      if(kernel.size() % 2 != 1) {
        BRICK_THROW(brick::common::ValueException, "filterRows()",
                    "Argument kernel must have an odd number of elements.");
      }

      switch(kernel.size()) {
      case 3:
        privateCode::filterRows_3(outputArray, inputArray, kernel);
        break;
      case 5:
        privateCode::filterRows_5(outputArray, inputArray, kernel);
        break;
      case 7:
        privateCode::filterRows_7(outputArray, inputArray, kernel);
        break;
      case 9:
        privateCode::filterRows_9(outputArray, inputArray, kernel);
        break;
      case 11:
        privateCode::filterRows_11(outputArray, inputArray, kernel);
        break;
      default:
        privateCode::filterRows_generic(outputArray, inputArray, kernel);
        break;
      }
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_FILTER_IMPL_HH */
