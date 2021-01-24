/**
***************************************************************************
* @file brick/utilities/pythonIO.hh
*
* Header file declaring routines for importing & exporting things in
* a format that's readable & writable by the python programming language.
*
* Copyright (C) 2005-2011, David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_PYTHONIO_HH
#define BRICK_UTILITIES_PYTHONIO_HH

#include <iostream>
#include <string>

namespace brick {

  namespace utilities {

    /**
     * This function formats a 1D array in a way that python can parse
     * into a built-in list.  The output will look something like this:
     *
     *   [1, 2, 3, 4]
     *
     * The array object must provide stl-style begin() and end()
     * member functions, and an ArrayType::const_iter typedef.
     *
     * @param array0 This argument is the array to be formatted.
     *
     * @return The return value is a string containing the
     * python-friendly text.
     */
    template<class ArrayType>
    std::string
    toPythonList(const ArrayType& array0);


    /**
     * This function formats a 1D Array object in a way that python can
     * parse into a Numpy array.  The output will look something like
     * this:
     *
     *   array([1, 2, 3, 4])
     *
     * The array object must provide stl-style begin() and end()
     * member functions.
     *
     * @param array0 This argument is the array to be formatted.
     *
     * @return The return value is a string containing the
     * python-friendly text.
     */
    template<class ArrayType>
    std::string
    toPythonNumeric1DArray(const ArrayType& array0);


    /**
     * This function formats an 2D Array object in a way that python can
     * parse into a Numeric array.  The output will look something like
     * this:
     *
     *   array([[1, 2], [3, 4]])
     *
     * The provided array object should provide the following member
     * functions.
     *
     * - rows() const: This member function should return the total
     *   number of rows in the array.
     *
     * - row(size_t) const: This member function should return an
     *   array-like object that addresses a specific row of the array.
     *   This object must be compatible with the function toPythonList(),
     *   which is declared above.
     *
     * @param array0 This argument is the array to be formatted.
     *
     * @return The return value is a string containing the
     * python-friendly text.
     */
    template<class ArrayType>
    std::string
    toPythonNumeric2DArray(const ArrayType& array0);

  } // namespace utilities

} // namespace brick


/* ================ Function definitions follow ================ */

namespace brick {

  namespace utilities {

    template<class ArrayType>
    std::string
    toPythonList(const ArrayType& array0)
    {
      // We'll use an ostringstream to format the string.
      std::ostringstream outputStream;

      // Python list syntax starts with a square-bracket.
      outputStream << "[";

      // We'll iterate over every element in the array.
      size_t count = 0;
      typename ArrayType::const_iter iter = array0.begin();
      while(iter != array0.end()) {

        // This variable will keep track of how many entries have been
        // written to each line, so that we can politely insert carriage
        // returns.
        size_t lineCount = 0;

        // Loop just enough to complete one line.
        while((iter != array0.end()) && (lineCount < 20)) {

          // If it's not the first entry, we need separating punctuation.
          if(count != 0) {
            // Usually a comma is all that's necessary.
            outputStream << ", ";

            // If we just finished a line, we'll need a carriage return and
            // an indent.
            if(lineCount == 0) {
              outputStream << "\n\t\t";
            }
          }

          // Write the next element.
          outputStream << *iter;

          // Update the counters.
          ++count;
          ++lineCount;
          ++iter;
        }
      }

      // Add the closing square-bracket.
      outputStream << "]";

      // And return the formatted string.
      return outputStream.str();
    }


    template<class ArrayType>
    std::string
    toPythonNumeric1DArray(const ArrayType& array0)
    {
      // We'll use an ostringstream to format the string.
      std::ostringstream outputStream;

      // Start with the appropriate syntax for a Numeric array,
      // and the the opening square-bracket of the top-level list.
      outputStream << "array(";
      outputStream << toPythonList(array0);

      // Close the open syntax.
      outputStream << ")\n";

      // And return the formatted string.
      return outputStream.str();
    }


    template<class ArrayType>
    std::string
    toPythonNumeric2DArray(const ArrayType& array0)
    {
      // We'll use an ostringstream to format the string.
      std::ostringstream outputStream;

      // Start with the appropriate syntax for a Numeric array,
      // and the the opening square-bracket of the top-level list.
      outputStream << "array([";

      // Now iterate over each row.
      for(size_t index0 = 0; index0 < array0.rows(); ++index0) {
        if(index0 != 0) {
          outputStream << ",\n\t";
        }
        outputStream << toPythonList(array0.row(index0));
      }

      // Close the open syntax.
      outputStream << "])\n";

      // And return the formatted string.
      return outputStream.str();
    }

  } // namespace utilities

} // namespace brick

#endif /* #ifndef BRICK_UTILITIES_PYTHONIO_HH */
