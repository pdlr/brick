/**
***************************************************************************
* @file brick/numeric/utilities.hh
*
* Header file declaring useful functions for dlrNumeric.
*
* Copyright (C) 2001-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_UTILITIES_IMPL_HH
#define BRICK_NUMERIC_UTILITIES_IMPL_HH

// This file is included by utilites.hh, and should not be directly
// included by user code, so no need to include utilities.hh here.
//
// #include <brick/numeric/utilites.hh>


#include <cmath>
#include <sstream>
#include <numeric>

#include <brick/numeric/mathFunctions.hh>
#include <brick/numeric/differentiableScalar.hh>

namespace brick {

  namespace numeric {

    /**
     ** @cond privateCode
     **
     ** "Private" namespace containing some functors to make
     ** implementation easier in utilities.h.  Note, since this
     ** (regrettably) has to be part of the .h file, it's important not
     ** to put anything in this namespace which is not either a template
     ** or an inline function.  If you violate this rule, you will
     ** probably get multiply defined symbols when you link your
     ** program.
     **/
    namespace privateCode {

      /**
       ** This is the generic absFunctor implementation.  The purpose of
       ** absFunctor is to compute the absolute value of the
       ** operator()(Type const&) input argument.
       **/
      template <class Type>
      struct absFunctor : public std::unary_function<Type, Type> {
        /**
         * This is the operator which does all the work.  Since there's
         * no generic way to compute the absolute value of a type, the
         * generic implementation simply throws an exception.
         *
         * @return This function throws an exception, and so does not return.
         */
        inline Type operator()(Type const& input) {
          BRICK_THROW(brick::common::NotImplementedException,
                     "absFunctor<Type>::operator()(Type const&)",
                     "absFunctor must be specialized for each type.");
          return static_cast<Type>(0);
        }
      };

      /**
       ** This functor specializes absFunctor for long double floats.
       **/
      template <>
      struct absFunctor<long double>
        : public std::unary_function<long double, long double> {
        /**
         * This is the operator which does all the work.  It returns the
         * absolute value of its input argument.
         *
         * @return This function throws an exception, and so does not return.
         */
        inline double operator()(long double const& input) {
          // fabsl doesn't appear to be commonly available.
          // return std::fabsl(input);
          return (input > 0.0) ? input : -input;
        }
      };


      /**
       ** This functor specializes absFunctor for double precision floats.
       **/
      template <>
      struct absFunctor<double> : public std::unary_function<double, double> {
        /**
         * This is the operator which does all the work.  It returns the
         * absolute value of its input argument.
         *
         * @return This function throws an exception, and so does not return.
         */
        inline double operator()(double const& input) {
          return std::fabs(input);
        }
      };


      /**
       ** This functor specializes absFunctor for single precision floats.
       **/
      template <>
      struct absFunctor<float> : public std::unary_function<float, float> {
        /**
         * This is the operator which does all the work.  It returns the
         * absolute value of its input argument.
         *
         * @return This function throws an exception, and so does not return.
         */
        inline float operator()(float const& input) {
          // Strange.  fabsf shows up in the global namespace.  Is this
          // a bug in the g++ 3.3 stand library?
          return fabsf(input);
        }
      };


      /**
       ** This functor specializes absFunctor for long ints.
       **/
      template <>
      struct absFunctor<long int>
        : public std::unary_function<long int, long int> {
        /**
         * This is the operator which does all the work.  It returns the
         * absolute value of its input argument.
         *
         * @return This function throws an exception, and so does not return.
         */
        inline long int operator()(long int const& input) {
          return std::labs(input);
        }
      };


      /**
       ** This functor specializes absFunctor for long long ints.
       **/
      template <>
      struct absFunctor<long long int>
        : public std::unary_function<long long int, long long int> {
        /**
         * This is the operator which does all the work.  It returns the
         * absolute value of its input argument.
         *
         * @return This function throws an exception, and so does not return.
         */
        inline long long int operator()(long long int const& input) {
          // Many compilers don't support C99.
          // return std::llabs(input);
          return input >= 0LL ? input : -input;
        }
      };


      /**
       ** This functor specializes absFunctor for ints.
       **/
      template <>
      struct absFunctor<int> : public std::unary_function<int, int> {
        /**
         * This is the operator which does all the work.  It returns the
         * absolute value of its input argument.
         *
         * @return This function throws an exception, and so does not return.
         */
        inline int operator()(int const& input) {
          return std::abs(input);
        }
      };


    }
    /** @endcond **/


    // This function returns an array of the same size and element type
    // as its input argument, in which each element is set to the
    // absolute value of the corresponding element of the input array.
    template <class Type>
    Array1D<Type>
    abs(Array1D<Type> const& array0)
    {
      Array1D<Type> result(array0.size());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     privateCode::absFunctor<Type>());
      return result;
    }


    // This function returns an array of the same size and element type
    // as its input argument, in which each element is set to the
    // absolute value of the corresponding element of the input array.
    template <class Type>
    Array2D<Type>
    abs(Array2D<Type> const& array0)
    {
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     privateCode::absFunctor<Type>());
      return result;
    }


    // This function returns true if each element of its argument is
    // false, and returns false otherwise.
    template <class Type>
    bool
    allFalse(Array1D<Type> const& array0)
    {
      for(typename Array1D<Type>::const_iterator iter = array0.begin();
          iter != array0.end();
          ++iter) {
        if(*iter) {
          return false;
        }
      }
      return true;
    }


    // This function returns true if each element of its argument is
    // true, and returns false otherwise.
    template <class Type>
    bool
    allTrue(Array1D<Type> const& array0)
    {
      for(typename Array1D<Type>::const_iterator iter = array0.begin();
          iter != array0.end();
          ++iter) {
        if(!(*iter)) {
          return false;
        }
      }
      return true;
    }


    // This function returns true if any element of its argument is
    // false, and returns false otherwise.
    template <class Type>
    inline bool
    anyFalse(Array1D<Type> const& array0)
    {
      return !allTrue(array0);
    }


    // This function returns true if any element of its argument is
    // true, and returns false otherwise.
    template <class Type>
    inline bool
    anyTrue(Array1D<Type> const& array0)
    {
      return !allFalse(array0);
    }


    // This function returns the index of the largest element of its input
    // array.
    template <class Type>
    inline size_t
    argmax(Array1D<Type> const& array0)
    {
      return argmax(array0, std::less<Type>());
    }


    // This function returns the index of the largest element of its input
    // sequence.
    template <class IterType>
    inline size_t
    argmax(IterType beginIter, IterType endIter)
    {
      return argmax(beginIter, endIter,
                    std::less<typename IterType::value_type>());
    }


    // This function returns the index of the largest element of its
    // input array, where largeness is defined by the second argument.
    template <class Type, class Functor>
    inline size_t
    argmax(Array1D<Type> const& array0, Functor comparator)
    {
      return argmax(array0.begin(), array0.end(), comparator);
    }


    // This function returns the index of the largest element of its
    // input sequence, where largeness is defined by the second argument.
    template <class IterType, class Functor>
    inline size_t
    argmax(IterType beginIter, IterType endIter, Functor comparator)
    {
      return static_cast<size_t>(
        std::max_element(beginIter, endIter, comparator) - beginIter);
    }


    // This function returns an Index2D instance indicating which is
    // the largest element of its input array.
    template <class Type>
    inline Index2D
    argmax2D(Array2D<Type> const& array0)
    {
      return argmax2D(array0, std::less<Type>());
    }


    // This function returns an Index2D instance indicating which is
    // the largest element of its input array, where largeness is
    // defined by the second argument.
    template <class Type, class Functor>
    inline Index2D
    argmax2D(Array2D<Type> const& array0, Functor comparator)
    {
      size_t ravelIndex = argmax(array0.ravel(), comparator);
      size_t row = ravelIndex / array0.columns(); // Int division.
      size_t column = ravelIndex - row * array0.columns();
      return Index2D(row, column);
    }


    // This function returns the index of the smallest element of its input
    // array.
    template <class Type>
    inline size_t
    argmin(Array1D<Type> const& array0)
    {
      return argmin(array0, std::less<Type>());
    }


    // This function returns the index of the smallest element of its
    // input array, where largeness is defined by the second argument.
    template <class Type, class Functor>
    inline size_t
    argmin(Array1D<Type> const& array0, Functor comparator)
    {
      Type const* minPtr = std::min_element(array0.begin(), array0.end(),
                                            comparator);
      return static_cast<size_t>(minPtr - array0.begin());
    }


    // This function returns an Index2D instance indicating which is
    // the smallest element of its input array.
    template <class Type>
    inline Index2D
    argmin2D(Array2D<Type> const& array0)
    {
      return argmin2D(array0, std::less<Type>());
    }


    // This function returns an Index2D instance indicating which is
    // the smallest element of its input array, where smallness is
    // defined by the second argument.
    template <class Type, class Functor>
    inline Index2D
    argmin2D(Array2D<Type> const& array0, Functor comparator)
    {
      size_t ravelIndex = argmin(array0.ravel(), comparator);
      size_t row = ravelIndex / array0.columns(); // Int division.
      size_t column = ravelIndex - row * array0.columns();
      return Index2D(row, column);
    }


    // This function returns an array of indices, result, so that the
    // sequence (array0[result[0]], array0[result[1]],
    // array0[result[2]], ...) is sorted from smallest to largest using
    // operator<().
    template <class Type>
    Array1D<size_t>
    argsort(Array1D<Type> const& array0)
    {
      Array1D< std::pair<Type, size_t> > sortVector(array0.size());
      Array1D<size_t> resultVector(array0.size());

      for(size_t index0 = 0; index0 < array0.size(); ++index0) {
        sortVector[index0] = std::pair<Type, size_t>(array0[index0], index0);
      }
      std::sort(sortVector.begin(), sortVector.end());
      std::transform(sortVector.begin(), sortVector.end(), resultVector.begin(),
                     ExtractSecondFunctor<Type, size_t>());
      return resultVector;
    }


//   // This function returns an array of indices, result, so that the
//   // sequence (array0[result[0]], array0[result[1]],
//   // array0[result[2]], ...) is sorted from smallest to largest using
//   // the supplied comparison operator.
//   template <class Type, class Functor>
//   Array1D<size_t>
//   argsort(Array1D<Type> const& array0, Functor comparator)
//   {
//     Array1D< std::pair<Type, size_t> > sortVector(array0.size());
//     Array1D<size_t> resultVector(array0.size());

//     for(size_t index0 = 0; index0 < array0.size(); ++index0) {
//       sortVector[index0] = std::pair<Type, size_t>(array0[index0], index0);
//     }
//     std::sort(sortVector.begin(), sortVector.end(), comparator);
//     std::transform(sortVector.begin(), sortVector.end(), resultVector.begin(),
//                    ExtractSecondFunctor<Type, size_t>());
//     return resultVector;
//   }


    // This function returns an Array1D in which each element has the
    // value of the largest element in one row or column of the input
    // Array2D.
    template <class Type>
    inline Array1D<Type>
    axisMaximum(Array2D<Type> const& array0, size_t axis)
    {
      return axisMaximum(array0, axis, std::less<Type>());
    }


    // This function returns an Array1D in which each element has the
    // value of the largest element in one row or column of the input
    // Array2D, where largeness is defined by the third argument.
    template <class Type, class Functor>
    Array1D<Type>
    axisMaximum(Array2D<Type> const& array0, size_t axis, Functor comparator)
    {
      Array1D<Type> result;
      switch(axis) {
      case 0:
        result.reinit(array0.columns());
        for(size_t column = 0; column < array0.columns(); ++column) {
          Type const* dataPtr = array0.data(0, column);
          Type columnMax = *dataPtr;
          size_t stride = array0.columns();
          for(size_t row = 0; row < array0.rows(); ++row) {
            if(!comparator(*dataPtr, columnMax)) {
              columnMax = *dataPtr;
            }
            dataPtr += stride;
          }
          result(column) = columnMax;
        }
        break;
      case 1:
        result.reinit(array0.rows());
        for(size_t row = 0; row < array0.rows(); ++row) {
          Type const* dataPtr = array0.data(row, 0);
          Type rowMax = *dataPtr;
          for(size_t column = 0; column < array0.columns(); ++column) {
            if(!comparator(*dataPtr, rowMax)) {
              rowMax = *dataPtr;
            }
            ++dataPtr;
          }
          result(row) = rowMax;
        }
        break;
      default:
        std::ostringstream message;
        message << "Axis " << axis << " is invalid for an Array2D.";
        BRICK_THROW(brick::common::IndexException, "axisMaximum(Array2D const&, size_t, ...)",
                   message.str().c_str());
        break;
      }
      return result;
    }


    // This function returns an Array1D in which each element has the
    // value of the smallest element in one row or column of the input
    // Array2D.
    template <class Type>
    inline Array1D<Type>
    axisMinimum(Array2D<Type> const& array0, size_t axis)
    {
      return axisMinimum(array0, axis, std::less<Type>());
    }


    // This function returns an Array1D in which each element has the
    // value of the smallest element in one row or column of the input
    // Array2D, where smallness is defined by the third argument.
    template <class Type, class Functor>
    Array1D<Type>
    axisMinimum(Array2D<Type> const& array0, size_t axis, Functor comparator)
    {
      Array1D<Type> result;
      switch(axis) {
      case 0:
        result.reinit(array0.columns());
        for(size_t column = 0; column < array0.columns(); ++column) {
          Type const* dataPtr = array0.data(0, column);
          Type columnMax = *dataPtr;
          size_t stride = array0.columns();
          for(size_t row = 0; row < array0.rows(); ++row) {
            if(comparator(*dataPtr, columnMax)) {
              columnMax = *dataPtr;
            }
            dataPtr += stride;
          }
          result(column) = columnMax;
        }
        break;
      case 1:
        result.reinit(array0.rows());
        for(size_t row = 0; row < array0.rows(); ++row) {
          Type const* dataPtr = array0.data(row, 0);
          Type rowMax = *dataPtr;
          for(size_t column = 0; column < array0.columns(); ++column) {
            if(comparator(*dataPtr, rowMax)) {
              rowMax = *dataPtr;
            }
            ++dataPtr;
          }
          result(row) = rowMax;
        }
        break;
      default:
        std::ostringstream message;
        message << "Axis " << axis << " is invalid for an Array2D.";
        BRICK_THROW(brick::common::IndexException,
                    "axisMinimum(Array2D const&, size_t, ...)",
                   message.str().c_str());
        break;
      }
      return result;
    }


    // This function returns an Array1D in which each element has the
    // sum of one row or column of the input Array2D.
    template <class ResultType, class Type>
    inline Array1D<ResultType>
    axisSum(Array2D<Type> const& array0, size_t axis)
    {
      return axisSum<ResultType>(
        array0, axis, static_cast<ResultType>(0), std::plus<ResultType>());
    }


    // This function returns an Array1D in which each element has the
    // sum of one row or column of the input Array2D.
    template <class ResultType, class Type, class Functor>
    Array1D<ResultType>
    axisSum(Array2D<Type> const& array0, size_t axis,
            ResultType const& initialValue, Functor adder)
    {
      Array1D<ResultType> result;
      switch(axis) {
      case 0:
        result.reinit(array0.columns());
        for(size_t column = 0; column < array0.columns(); ++column) {
          ResultType columnSum = initialValue;
          Type const* dataPtr = array0.data(0, column);
          size_t stride = array0.columns();
          for(size_t row = 0; row < array0.rows(); ++row) {
            columnSum = adder(columnSum, *dataPtr);
            dataPtr += stride;
          }
          result(column) = columnSum;
        }
        break;
      case 1:
        result.reinit(array0.rows());
        for(size_t row = 0; row < array0.rows(); ++row) {
          ResultType rowSum = initialValue;
          Type const* dataPtr = array0.data(row, 0);
          for(size_t column = 0; column < array0.columns(); ++column) {
            rowSum = adder(rowSum, *dataPtr);
            ++dataPtr;
          }
          result(row) = rowSum;
        }
        break;
      default:
        std::ostringstream message;
        message << "Axis " << axis << " is invalid for an Array2D.";
        BRICK_THROW(brick::common::IndexException, "axisSum(Array2D const&, size_t, ...)",
                   message.str().c_str());
        break;
      }
      return result;
    }


    // This function swaps values between its two arguments so that
    // the value of each element of corner0 is less than or equal to
    // the corresponding value of corner1.
    template <class Type>
    bool
    cleanupCorners(Vector2D<Type>& corner0, Vector2D<Type>& corner1)
    {
      bool returnValue = true;

      if(corner0.getX() > corner1.getX()) {
        Type maxX = corner0.getX();
        corner0.setValue(corner1.getX(), corner0.getY());
        corner1.setValue(maxX, corner1.getY());
      } else if(corner0.getX() == corner1.getX()) {
        returnValue = false;
      }

      if(corner0.getY() > corner1.getY()) {
        Type maxY = corner0.getY();
        corner0.setValue(corner0.getX(), corner1.getY());
        corner1.setValue(corner1.getX(), maxY);
      } else if(corner0.getY() == corner1.getY()) {
        returnValue = false;
      }

      return returnValue;
    }


    // columnIndices(rows, columns): Returns an Array2D in which each
    // element contains the index of its column.
    template <class Type>
    Array2D<Type>
    columnIndices(size_t rows, size_t columns)
    {
      Array2D<Type> returnValue(rows, columns);
      Array1D<Type> modelRow = range(static_cast<Type>(0),
                                     static_cast<Type>(columns),
                                     static_cast<Type>(1));
      for(size_t row=0; row < rows; ++row) {
        std::copy(modelRow.begin(), modelRow.end(), returnValue.data(row, 0));
      }
      return returnValue;
    }


    // This function selects those elements of an input Array1D which
    // correspond to "true" values of a mask Array1D, and returns an
    // Array1D containing only those elements.
    template <class Type0, class Type1>
    inline Array1D<Type1>
    compress(Array1D<Type0> const& condition,
             Array1D<Type1> const& input)
    {
      size_t numTrue = count(condition);
      return compress(condition, input, numTrue);
    }


    // This function behaves in exactly the same way as compress(const
    // Array1D&, Array1D const&), above, but it permits the user to
    // specify the number of true elements in the condition array.
    template <class Type0, class Type1>
    Array1D<Type1>
    compress(Array1D<Type0> const& condition,
             Array1D<Type1> const& input,
             size_t numTrue)
    {
      if(condition.size() != input.size()) {
        std::ostringstream message;
        message << "Condition and input arguments must have the same "
                << "size, but condition has size = " << condition.size()
                << ", while input has size = " << input.size() << ".";
        BRICK_THROW(brick::common::ValueException,
                   "compress(Array1D const&, Array1D const&, size_t)",
                   message.str().c_str());
      }
      Array1D<Type1> result(numTrue);
      // copy_mask(input.begin(), input.end(), condition.begin(),
      //           result.begin());
      //
      // Hmm... copy_mask was previously defined in dlrLib3 as follows:
      //
      // template <class In0, class In1, class Out>
      // Out copy_mask(In0 first, In0 last, In1 mask, Out res)
      // {
      //   while(first != last) {
      //     if(*mask) {
      //       *res++ = *first;
      //     }
      //     ++mask;
      //     ++first;
      //   }
      //   return res;
      // }

      typename Array1D<Type1>::const_iterator first = input.begin();
      typename Array1D<Type1>::const_iterator last = input.end();
      typename Array1D<Type0>::const_iterator mask = condition.begin();
      typename Array1D<Type1>::iterator target = result.begin();
      while(first != last) {
        if(*mask) {
          *(target++) = *first;
        }
        ++mask;
        ++first;
      }

      return result;
    }

    // This element counts the number of elements of the input array
    // which evaluate to true, and returns that number.
    template <class Type>
    inline size_t
    count(Array1D<Type> const& x)
    {
      return std::count_if(x.begin(), x.end(), StaticCastFunctor<Type, bool>());
    }

    // This function computes the cross product of two Vector3D
    // instances.
    template <class Type>
    inline Vector3D<Type>
    cross(Vector3D<Type> const& vector0, Vector3D<Type> const& vector1) {
      return Vector3D<Type>(
        (vector0.y() * vector1.z()) - (vector0.z() * vector1.y()),
        (vector0.z() * vector1.x()) - (vector0.x() * vector1.z()),
        (vector0.x() * vector1.y()) - (vector0.y() * vector1.x()));
    }


    // This function computes the inner product of two input arrays, and
    // allows the user to control which type is used to do the
    // calculation.
    template <class Type2, class Type1, class Type0>
    inline Type2
    dot(Array1D<Type0> const& array0, Array1D<Type1> const& array1)
    {
      if(array0.size() != array1.size()) {
        std::ostringstream message;
        message << "Input arguments must have matching sizes, "
                << "but array0.size() == " << array0.size()
                << " and array1.size() == " << array1.size() << "."
                << std::endl;
        BRICK_THROW(brick::common::ValueException, "dot()", message.str().c_str());
      }
      return std::inner_product(array0.begin(), array0.end(), array1.begin(),
                                static_cast<Type2>(0));
    }


    // This function computes the inner product of two Vector2D instances.
    template <class Type1, class Type0>
    inline Type1
    dot(Vector2D<Type0> const& vector0, Vector2D<Type0> const& vector1)
    {
      return ((static_cast<Type1>(vector0.x())
               * static_cast<Type1>(vector1.x()))
              + (static_cast<Type1>(vector0.y())
                 * static_cast<Type1>(vector1.y())));
    }


    // This function computes the inner product of two Vector3D instances.
    template <class Type1, class Type0>
    inline Type1
    dot(Vector3D<Type0> const& vector0, Vector3D<Type0> const& vector1)
    {
      return ((static_cast<Type1>(vector0.x())
               * static_cast<Type1>(vector1.x()))
              + (static_cast<Type1>(vector0.y())
                 * static_cast<Type1>(vector1.y()))
              + (static_cast<Type1>(vector0.z())
                 * static_cast<Type1>(vector1.z())));
    }


    // This function computes the matrix X such that A * x = X * vec(A).
    template <class Type>
    Array2D<Type>
    equivalentMatrix(Array1D<Type> const& vector0, size_t rowsInMatrix)
    {
      Array2D<Type> XMatrix = zeros<Type>(
        rowsInMatrix, vector0.size() * rowsInMatrix);
      for(size_t XRow = 0; XRow != rowsInMatrix; ++XRow) {
        for(size_t index0 = 0; index0 < vector0.length(); ++index0) {
          XMatrix(XRow, index0 * rowsInMatrix + XRow) = vector0(index0);
        }
      }
      return XMatrix;
    }

    // This function returns a 2D matrix with the specified shape in
    // which the elements on the diagonal are set to 1, and all other
    // elements are set to 0.
    template <class Type>
    Array2D<Type>
    identity(size_t rows, size_t columns)
    {
      Array2D<Type> returnArray = zeros<Type>(rows, columns);
      size_t rank = std::min(rows, columns);
      for(size_t index = 0; index < rank; ++index) {
        returnArray(index, index) = static_cast<Type>(1);
      }
      return returnArray;
    }


    // This function returns an array in which each element is the
    // natural logarithm of the corresponding element of its input.
    template <class Type>
    Array1D<Type>
    ln(Array1D<Type> const& array0)
    {
      Array1D<Type> result(array0.size());
      // std::transform(array0.begin(), array0.end(), std::log);
      for(size_t index0 = 0; index0 < array0.size(); ++index0) {
        result[index0] = std::log(array0[index0]);
      }
      return result;
    }


//   // Strange troubles compiling these next two functions.
//
//   // This function returns an array in which each element is the
//   // natural logarithm of the corresponding element of its input.
//   template <>
//   Array1D<float>
//   ln(Array1D<float> const& array0)
//   {
//     Array1D<float> result(array0.size());
//     std::transform(array0.begin(), array0.end(), std::logf);
//     return result;
//   }


//   // This function returns an array in which each element is the
//   // natural logarithm of the corresponding element of its input.
//   template <>
//   Array1D<long double>
//   ln(Array1D<long double> const& array0)
//   {
//     Array1D<long double> result(array0.size());
//     std::transform(array0.begin(), array0.end(), std::logl);
//     return result;
//   }


    // This function returns an array in which each element is the
    // natural logarithm of the corresponding element of its input.
    template <class Type>
    Array2D<Type>
    ln(Array2D<Type> const& array0)
    {
      Array2D<Type> result(array0.size());
      // std::transform(array0.begin(), array0.end(), std::log);
      for(size_t index0 = 0; index0 < array0.size(); ++index0) {
        result[index0] = std::log(array0[index0]);
      }
      return result;
    }


    // This function returns an Array2D instance in which the value of
    // each element is the logical not of the corresponding element of
    // the input array.
    template <class Type>
    Array2D<Type>
    logicalNot(Array2D<Type> const& array0)
    {
      Array2D<Type> result(array0.rows(), array0.columns());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     unaryComposeFunctor(StaticCastFunctor<bool, Type>(),
                                         std::logical_not<Type>()));
      return result;
    }


    // This function computes the magnitude of its input argument.
    template <class Type1, class Type0>
    inline Type1
    magnitude(Array1D<Type0> const& array0)
    {
      // TBD.
      return brick::numeric::squareRoot(magnitudeSquared<Type1>(array0));
    }


    // This function computes the magnitude of its input argument.
    template <class Type1, class Type0>
    inline Type1
    magnitude(Vector2D<Type0> const& vector0) {
      // TBD.
      return brick::numeric::squareRoot(magnitudeSquared<Type1>(vector0));
    }


    // This function computes the magnitude of its input argument.
    template <class Type1, class Type0>
    inline Type1
    magnitude(Vector3D<Type0> const& vector0) {
      return brick::numeric::squareRoot(magnitudeSquared<Type1>(vector0));
    }


    // This function computes the square of the magnitude of its input
    // argument.
    template <class Type1, class Type0>
    inline Type1
    magnitudeSquared(Array1D<Type0> const& array0)
    {
      return std::inner_product(array0.begin(), array0.end(),
                                array0.begin(), static_cast<Type0>(0));
    }


    // This function computes the square of the magnitude of its input
    // argument.
    template <class Type1, class Type0>
    inline Type1
    magnitudeSquared(Vector2D<Type0> const& vector0) {
      return ((static_cast<Type1>(vector0.x())
               * static_cast<Type1>(vector0.x()))
              + (static_cast<Type1>(vector0.y())
                 * static_cast<Type1>(vector0.y())));
    }


    // This function computes the square of the magnitude of its input
    // argument.
    template <class Type1, class Type0>
    inline Type1
    magnitudeSquared(Vector3D<Type0> const& vector0) {
      return ((static_cast<Type1>(vector0.x())
               * static_cast<Type1>(vector0.x()))
              + (static_cast<Type1>(vector0.y())
                 * static_cast<Type1>(vector0.y()))
              + (static_cast<Type1>(vector0.z())
                 * static_cast<Type1>(vector0.z())));
    }


    // This function computes a vector * matrix product just like
    // matrixMultiply(Array1D<Type> const&, Array2D<Type> const&), above.
    // This function differs in that the element type of the return
    // value is set explicitly using a third argument.
    template <class Type2, class Type1, class Type0>
    Array1D<Type2>
    matrixMultiply(Array1D<Type0> const& vector0, Array2D<Type1> const& matrix0)
    {
      if(vector0.size() != matrix0.rows()) {
        std::ostringstream message;
        message << "Can't left-multiply a "
                << matrix0.rows() << " x " << matrix0.columns()
                << " matrix by a " << vector0.size() << " element vector\n";
        BRICK_THROW(brick::common::ValueException, "matrixMultiply()",
                    message.str().c_str());
      }
      Array1D<Type2> result = zeros<Type2>(matrix0.columns());
      size_t index = 0;
      for(size_t row = 0; row < matrix0.rows(); ++row) {
        for(size_t column = 0; column < matrix0.columns(); ++column) {
          result[column] += (static_cast<Type2>(vector0[row])
                             * static_cast<Type2>(matrix0[index]));
          ++index;
        }
      }
      return result;
    }


    // This function computes a matrix * vector product just like
    // matrixMultiply(Array2D<Type> const&, Array1D<Type> const&), above.
    // This function differs in that the element type of the return
    // value is set explicitly using a third argument.
    template <class Type2, class Type1, class Type0>
    Array1D<Type2>
    matrixMultiply(Array2D<Type0> const& matrix0, Array1D<Type1> const& vector0)
    {
      if(vector0.size() != matrix0.columns()) {
        std::ostringstream message;
        message << "matrixMultiply() -- can't right-multiply a "
                << matrix0.rows() << " x " << matrix0.columns()
                << " matrix by a " << vector0.size() << " element vector\n";
        BRICK_THROW(brick::common::ValueException, "matrixMultiply()", message.str().c_str());
      }
      Array1D<Type2> result(matrix0.rows());
      for(size_t row = 0; row < matrix0.rows(); ++row) {
        result[row] = dot<Type2>(vector0, matrix0.row(row));
      }
      return result;
    }


    // This function computes a matrix * matrix product, just like
    // matrixMultiply(Array2D<Type> const&, Array2D<Type> const&), above.
    // This function differs in that the element type of the return
    // value is set explicitly using a third argument.
    template <class Type2, class Type1, class Type0>
    Array2D<Type2>
    matrixMultiply(Array2D<Type0> const& matrix0, Array2D<Type1> const& matrix1)
    {
      if(matrix1.rows() != matrix0.columns()) {
        std::ostringstream message;
        message << "matrixMultiply() -- can't left-multiply a "
                << matrix1.rows() << " x " << matrix1.columns()
                << " matrix by a "
                << matrix0.rows() << " x " << matrix0.columns()
                << " matrix\n";
        BRICK_THROW(brick::common::ValueException, "matrixMultiply()", message.str().c_str());
      }
      Array2D<Type2> result = zeros<Type2>(matrix0.rows(), matrix1.columns());
      for(size_t resultRow = 0; resultRow < result.rows(); ++resultRow) {
        for(size_t resultColumn = 0; resultColumn < result.columns();
            ++resultColumn) {
          for(size_t index = 0; index < matrix0.columns(); ++index) {
            result(resultRow, resultColumn) +=
              (static_cast<Type2>(matrix0(resultRow, index))
               * static_cast<Type2>(matrix1(index, resultColumn)));
          }
        }
      }
      return result;
    }


    // This function returns a copy of the largest element in the input
    // Array1D instance.
    template <class Type>
    inline Type
    maximum(Array1D<Type> const& array0)
    {
      return maximum(array0, std::less<Type>());
    }


    // This function returns a copy of the largest element in the input
    // Array1D instance, where largeness is defined by the return value
    // of the second argument.
    template <class Type, class Functor>
    Type
    maximum(Array1D<Type> const& array0, Functor comparator)
    {
      if(array0.size() == 0) {
        BRICK_THROW(brick::common::ValueException, "maximum()",
                   "Can't find the maximum element of an empty array.");
      }
      return *std::max_element(array0.begin(), array0.end(), comparator);
    }


    // This function returns a copy of the largest element in the input
    // Array2D instance.
    template <class Type>
    inline Type
    maximum(Array2D<Type> const& array0)
    {
      return maximum(array0, std::less<Type>());
    }


    // This function returns a copy of the largest element in the input
    // Array2D instance, where largeness is defined by the return value
    // of the second argument.
    template <class Type, class Functor>
    Type
    maximum(Array2D<Type> const& array0, Functor comparator)
    {
      Type maxElement = 0;

      if(array0.size() == 0) {
        BRICK_THROW(brick::common::ValueException, "maximum()",
                   "Can't find the maximum element of an empty array.");
      }

      if(array0.isContiguous()) {
        maxElement = *std::max_element(array0.begin(), array0.end(),
                                       comparator);
      } else {
        for(unsigned int rr = 0; rr < array0.rows(); ++rr) {
          Type candidate = maximum(array0.getRow(rr), comparator);
          if(comparator(maxElement, candidate)) {
            maxElement = candidate;
          }
        }
      }
      return maxElement;
    }


    // This function computes the average value, or geometric mean, of
    // the elements of input sequence.
    template <class Type, class Iterator>
    Type
    mean(Iterator beginIter, Iterator endIter)
    {
      Type meanValue = 0;
      size_t count = 0;
      while(beginIter != endIter) {
        meanValue += *beginIter;
        ++count;
        ++beginIter;
      }
      if(count == 0) {
        count = 1;
      }
      return meanValue / count;
    }


    // This function computes the average value, or geometric mean, of
    // the elements of its argument, and allows the user to specify the
    // precision with which the computation is carried out.
    template <class Type1, class Type0>
    inline Type1
    mean(Array1D<Type0> const& array0)
    {
      return sum<Type1>(array0) / static_cast<Type1>(array0.size());
    }


    // This function returns the centroid of a 1D array.
    template <class FloatType, class Type>
    FloatType
    getCentroid(Array1D<Type> const& signal) {
      return getCentroid<FloatType>(signal.begin(), signal.end());
    }


    // This function returns the centroid of a 1D array.
    template <class FloatType, class IterType>
    FloatType
    getCentroid(IterType beginIter, IterType endIter)
    {
      FloatType weightedAccumulator = FloatType(0);
      FloatType accumulator = FloatType(0);
      std::size_t ii = 1;
      while(beginIter != endIter) {
        accumulator += *beginIter;
        weightedAccumulator += static_cast<FloatType>(ii) * *beginIter;
        ++ii;
        ++beginIter;
      }

      if(accumulator < NumericTraits<FloatType>::epsilon()) {
        BRICK_THROW(brick::common::ValueException,
                    "getCentroid()",
                    "Sum of input sequence is too small to reliably "
                    "compute centroid.");
      }

      // Remember that indices here are zero-based.
      return (weightedAccumulator / accumulator) - 1;
    }


    // This function estimates the mean and covariance of a set of
    // vectors, which are represented by the rows (or columns) of the
    // input 2D array.
    template <class Type>
    void
    getMeanAndCovariance(Array2D<Type> const& sampleArray,
                         Array1D<double>& meanArray,
                         Array2D<double>& covarianceArray,
                         size_t majorAxis)
    {
      // Check input.
      if(sampleArray.size() == 0) {
        meanArray.reinit(0);
        covarianceArray.reinit(0, 0);
        return;
      }

      // Compute number of samples and sample dimensionality.
      size_t dimensionality;
      size_t numberOfSamples;
      if(majorAxis == 0) {
        dimensionality = sampleArray.columns();
        numberOfSamples = sampleArray.rows();
      } else {
        dimensionality = sampleArray.rows();
        numberOfSamples = sampleArray.columns();
      }

      // Now compute the mean value of the sample points.
      // For efficiency, we simultaneously compute covariance.
      // We make use of the identity:
      //
      //   covariance = E[(x - u)*(x - u)^T] = E[x * x^T] - u * u^T
      //
      // where E[.] denotes expected value and u is the mean.
      //
      // Note that the sample mean is an unbiased estimator for u, but
      // substituting the sample mean for u in the covariance equation
      // gives a biased estimator of covariance.  We resolve this by
      // using the unbiased estimator:
      //
      //   covariance = (1/(N-1))sum(x * x^T) - (N/(N-1)) * uHat * uHat^T
      //
      // where uHat is the sample mean.

      meanArray.reinit(dimensionality);
      meanArray = 0.0;
      covarianceArray.reinit(dimensionality, dimensionality);
      covarianceArray = 0.0;
      typename Array2D<Type>::const_iterator sampleIterator =
        sampleArray.begin();

      if(majorAxis == 0) {
        for(size_t rowIndex = 0; rowIndex < sampleArray.rows(); ++rowIndex) {
          for(size_t columnIndex = 0; columnIndex < sampleArray.columns();
              ++columnIndex) {
            // Update accumulator for mean.
            meanArray[columnIndex] += *sampleIterator;

            // Update accumulator for covariance.  We only compute the top
            // half of the covariance matrix for now, since it's symmetric.
            typename Array2D<Type>::const_iterator sampleIterator2 =
              sampleIterator;
            for(size_t covarianceIndex = columnIndex;
                covarianceIndex < dimensionality;
                ++covarianceIndex) {
              // Note(xxx): Should fix this to run faster.
              covarianceArray(columnIndex, covarianceIndex) +=
                *sampleIterator * *sampleIterator2;
              ++sampleIterator2;
            }
            ++sampleIterator;
          }
        }
      } else {
        for(size_t rowIndex = 0; rowIndex < sampleArray.rows(); ++rowIndex) {
          for(size_t columnIndex = 0; columnIndex < sampleArray.columns();
              ++columnIndex) {
            // Update accumulator for mean.
            meanArray[rowIndex] += *sampleIterator;

            // Update accumulator for covariance.  We only compute the top
            // half of the covariance matrix for now, since it's symmetric.
            typename Array2D<Type>::const_iterator sampleIterator2 =
              sampleIterator;
            for(size_t covarianceIndex = rowIndex;
                covarianceIndex < dimensionality;
                ++covarianceIndex) {
              // Note(xxx): Should fix this to run faster.
              covarianceArray(rowIndex, covarianceIndex) +=
                *sampleIterator * *sampleIterator2;
              sampleIterator2 += numberOfSamples;
            }
            ++sampleIterator;
          }
        }
      }

      // Finish up the computation of the mean.
      meanArray /= static_cast<double>(numberOfSamples);

      // Finish up the computation of the covariance matrix.
      for(size_t index0 = 0; index0 < dimensionality; ++index0) {
        for(size_t index1 = index0; index1 < dimensionality; ++index1) {
          // Add correction to the top half of the covariance matrix.
          covarianceArray(index0, index1) -=
            numberOfSamples * meanArray[index0] * meanArray[index1];
          covarianceArray(index0, index1) /= (numberOfSamples - 1);

          // And copy to the bottom half.
          if(index0 != index1) {
            covarianceArray(index1, index0) = covarianceArray(index0, index1);
          }
        }
      }
    }


    // This function estimates the mean and variance of a sequence of
    // scalar values.
    template <class Iter, class Type>
    void
    getMeanAndVariance(Iter beginIter, Iter endIter,
                       Type& mean, Type& variance)
    {
      Type sumOfElements = 0;
      Type sumOfSquaredElements = 0;
      unsigned int numberOfElements = 0;

      while(beginIter != endIter) {
        Type elementValue = *beginIter;
        sumOfElements += elementValue;
        sumOfSquaredElements += elementValue * elementValue;
        ++numberOfElements;
        ++beginIter;
      }

      // Here we use the unbiased estimator for variance:
      //
      //   var = sum_i((x_i - mu)^2) / (N - 1),
      //
      // where sum_i() is summation over i, x_i is the i^th element of
      // the sequence, ^2 denotes raising to the power of two, and mu
      // is the estimated mean of the sequence.  Multiplying (x_i -
      // mu)^2 out, and separating terms, we get
      //
      //   var = ((1 / (N - 1)) * sum_i(x_i * x_i)
      //         - 2 * (mu / (N - 1)) * sum_i(x_i)
      //         + (N / (N - 1)) * mu^2.
      //
      // Recognizing that sum_i(x_i) is N * mu, and combining terms,
      // we get
      //
      //   var = (1 / (N - 1)) * sum_i(x_i * x_i)
      //         - (N / (N - 1)) * mu * mu,
      //
      //   var = (sum_i(x_i * x_i) - N * mu * mu) / (N - 1)
      //
      // which is implemented below, except that we replace N * mu
      // with sum_i(x_i), because we've already computed that.
      if(numberOfElements >= 2) {
        mean = sumOfElements / numberOfElements;
        variance = ((sumOfSquaredElements - sumOfElements * mean)
                    / (numberOfElements - 1));
      } else if(numberOfElements >= 1) {
        mean = sumOfElements /* / numberOfElements */;
        variance = Type(0.0);
      } else {
        mean = Type(0.0);
        variance = Type(0.0);
      }
    }

    // This function estimates the mean and variance of a sequence of
    // scalar values, discarding values that appear to be outliers.
    template <class Iter, class Type>
    void
    getMeanAndVarianceRobust(Iter beginIter, Iter endIter,
                             double requiredConfidence,
                             double inlierProportion,
                             Type& mean, Type& variance)
    {
      // Check arguments.
      if((requiredConfidence < 0.0) || (requiredConfidence >= 1.0)) {
        BRICK_THROW(brick::common::ValueException,
                  "getMeanAndVarianceRobust()",
                  "Probability value requiredConfidence is out of range.");
      }
      if(inlierProportion > 1.0) {inlierProportion = 1.0;}
      if(inlierProportion <= 0.0) {
        BRICK_THROW(brick::common::ValueException, "getMeanAndVarianceRobust()",
                    "Inlier proportion must be greater than 0.");
      }
      // Warning(xxx): Not implemented.
      BRICK_THROW(brick::common::NotImplementedException,
                  "getMeanAndVarianceRobust()",
                  "This function still needs to be written!");
    }


    // This function returns a copy of the smallest element in the input
    // Array1D instance.
    template <class Type>
    inline Type
    minimum(Array1D<Type> const& array0)
    {
      if(array0.size() == 0) {
        BRICK_THROW(brick::common::ValueException, "minimum()",
                   "Can't find the minimum element of an empty array.");
      }
      return minimum(array0, std::less<Type>());
    }


    // This function returns a copy of the smallest element in the input
    // Array1D instance, where largeness is defined by the value of the
    // second argument.
    template <class Type, class Functor>
    Type
    minimum(Array1D<Type> const& array0, Functor comparator)
    {
      return *std::min_element(array0.begin(), array0.end(), comparator);
    }


    // This function uses the iterative method of Newton and Raphson to
    // search for a zero crossing of the supplied functor.
    template <class FUNCTOR>
    double newtonRaphson(double startPoint, FUNCTOR objectiveFunction,
                         double epsilon, size_t maxIterations)
    {
      double argument = startPoint;
      for(size_t count = 0; count < maxIterations; ++count) {
        double result = objectiveFunction(argument);
        if(std::fabs(result) < epsilon) {
          return argument;
        }
        double resultPrime = objectiveFunction.derivative(argument);
        if(resultPrime == 0.0) {
          BRICK_THROW(brick::common::ValueException, "newtonRaphson()", "Zero derivative.");
        }
        argument = argument - (result / resultPrime);
      }
      BRICK_THROW(brick::common::ValueException, "newtonRaphson()",
                "Root finding failed to converge.");
      return 0.0;  // keep compiler happy
    }


    // This function computes the normalized correlation of two Array1D
    // arguments.
    template <class Type2, class Type>
    Type2
    normalizedCorrelation(Array1D<Type> const& signal0,
                          Array1D<Type> const& signal1)
    {
      if(signal0.size() != signal1.size()) {
        BRICK_THROW(brick::common::ValueException, "normalizedCorrelation()",
                  "Input arrays must have the same size.");
      }
      Type2 oneOverN =
        static_cast<Type2>(1) / static_cast<Type2>(signal0.size());
      Type2 sum0 = static_cast<Type2>(0);
      Type2 sum1 = static_cast<Type2>(0);
      Type2 sum00 = static_cast<Type2>(0);
      Type2 sum01 = static_cast<Type2>(0);
      Type2 sum11 = static_cast<Type2>(0);
      typedef typename Array1D<Type>::const_iterator SignalIter;
      SignalIter iter0 = signal0.begin();
      SignalIter iter1 = signal1.begin();
      SignalIter end0 = signal0.end();
      while(iter0 != end0) {
        sum0 += *iter0;
        sum1 += *iter1;
        sum00 += (*iter0) * (*iter0);
        sum01 += (*iter0) * (*iter1);
        sum11 += (*iter1) * (*iter1);
        ++iter0;
        ++iter1;
      }
      Type2 numerator = sum01 - oneOverN * sum0 * sum1;
      Type2 denominator = brick::numeric::squareRoot(
        (sum00 - oneOverN * sum0 * sum0)
        * (sum11 - oneOverN * sum1 * sum1));
      return numerator / denominator;
    }


    // This function returns an Array1D of the specified size and type
    // in which the value of every element is initialized to 1.
    template <class Type>
    inline Array1D<Type>
    ones(int size)
    {
      Array1D<Type> result(size);
      result = static_cast<Type>(1);
      return result;
    }


    // This function returns an Array2D of the specified size and type
    // in which the value of every element is initialized to one.
    template <class Type>
    inline Array2D<Type>
    ones(int rows, int columns)
    {
      Array2D<Type> result(rows, columns);
      result = static_cast<Type>(1);
      return result;
    }


    // This function computes the outer product of two input Array1D
    // instances and allows the user to control which type is used to do the
    // computation.
    template <class Type2, class Type0, class Type1>
    Array2D<Type2>
    outerProduct(Array1D<Type0> const& x, Array1D<Type1> const& y)
    {
      Array2D<Type2> result(x.size(), y.size());
      typename Array2D<Type2>::iterator runningIterator = result.begin();
      typename Array1D<Type1>::const_iterator yBegin = y.begin();
      typename Array1D<Type1>::const_iterator yEnd = y.end();
      for(size_t index = 0; index < x.size(); ++index) {
        std::transform(yBegin, yEnd, runningIterator,
                       std::bind2nd(std::multiplies<Type2>(), x[index]));
        runningIterator += y.size();
      }
      return result;
    }


    // This function returns an Array1D in which the first element has
    // value equal to argument "start," and each subsequent element has
    // value equal to the previous element plus argument "stride."
    template <class Type>
    Array1D<Type>
    range(Type start, Type stop, Type stride)
    {
      if(stride == 0) {
        BRICK_THROW(brick::common::ValueException, "range(Type, Type, Type)",
                   "Argument \"stride\" must not be equal to 0.");
      }
      int length = static_cast<int>((stop - start) / stride);
      if(stride > 0) {
        if((length * stride) < (stop - start)) {
          ++length;
        }
      } else {
        if((length * stride) > (stop - start)) {
          ++length;
        }
      }
      Array1D<Type> result(length);
      for(int index = 0; index < length; ++index) {
        result(index) = start;
        start += stride;
      }
      return result;
    }


    // This function computes the RMS (Root Mean Square) value of the
    // elements of its argument.
    template <class Type1, class Type0>
    Type1
    rms(Array1D<Type0> const& array0)
    {
      Type1 accumulator = 0;
      for(size_t index = 0; index < array0.size(); ++index) {
        Type1 element = static_cast<Type1>(array0[index]);
        accumulator += element * element;
      }
      return brick::numeric::squareRoot(
        accumulator / static_cast<Type1>(array0.size()));
    }


    // rowIndices(rows, columns): Returns an Array2D in which each
    // element contains the index of its row.
    template <class Type>
    Array2D<Type>
    rowIndices(size_t rows, size_t columns)
    {
      Array2D<Type> returnValue(rows, columns);
      for(size_t row=0; row < rows; ++row) {
        Type* rowPtr = returnValue.data(row, 0);
        std::fill(rowPtr, rowPtr + columns, static_cast<Type>(row));
      }
      return returnValue;
    }


    // This function returns true if the two arrays have the same shape,
    // false otherwise.
    template <class Type0, class Type1>
    inline bool
    shapeMatch(Array1D<Type0> const& array0, Array1D<Type1> const& array1)
    {
      return array0.size() == array1.size();
    }


    // This function returns true if the two arrays have the same shape,
    // false otherwise.
    template <class Type0, class Type1>
    inline bool
    shapeMatch(Array2D<Type0> const& array0, Array2D<Type1> const& array1)
    {
      return ((array0.rows() == array1.rows())
              && (array0.columns() == array1.columns()));
    }


    // This function returns true if the two arrays have the same shape,
    // false otherwise.
    template <class Type0, class Type1>
    inline bool
    shapeMatch(Array3D<Type0> const& array0, Array3D<Type1> const& array1)
    {
      return ((array0.shape0() == array1.shape0())
              && (array0.shape1() == array1.shape1())
              && (array0.shape2() == array1.shape2()));
    }


    // skewSymmetric(x): Returns a skew symmetric matrix X such
    // that matrixMultiply(X, y) = cross(x, y).
    template <class Type>
    inline Array2D<Type>
    skewSymmetric(Array1D<Type> const& vector0)
    {
      if(vector0.size() != 3) {
        std::ostringstream message;
        message << "Argument must have exactly 3 elements, but instead has "
                << vector0.size() << "elements.";
        BRICK_THROW(brick::common::ValueException, "skewSymmetric()",
                    message.str().c_str());
      }
      Array2D<Type> returnVal(3, 3);
      returnVal(0, 0) = 0;
      returnVal(0, 1) = -vector0(2);
      returnVal(0, 2) = vector0(1);

      returnVal(1, 0) = vector0(2);
      returnVal(1, 1) = 0;
      returnVal(1, 2) = -vector0(0);

      returnVal(2, 0) = -vector0(1);
      returnVal(2, 1) = vector0(0);
      returnVal(2, 2) = 0;
      return returnVal;
    }


#if 0
    // This function computes the real roots of the quadratic polynomial
    // c0*x^2 + c1*x + c0 = 0.
    template <class Type>
    void
    solveQuadratic(Type c0, Type c1, Type c2,
                   Type& root0, Type& root1, bool& valid, bool)
    {
      valid = solveQuadratic(c0, c1, c2, root0, root1);
    }
#endif

    // This function computes the standard deviation of a group of
    // scalar samples.
    template <class Type0, class Type1>
    inline Type1
    standardDeviation(Array1D<Type0> const& array0)
    {
      return brick::numeric::squareRoot(variance<Type1>(array0));
    }


    // This function computes the sum of all elements in the input
    // array.
    template <class Type2, class Type>
    Type2
    sum(Array1D<Type> const& array0)
    {
      return std::accumulate(array0.begin(), array0.end(),
                             static_cast<Type2>(0));
    }


    // This function computes the sum of those elements of its
    // argument which lie within a rectangular region of interest.
    template <class Type2, class Type>
    Type2
    sum(Array2D<Type> const& array0,
        Index2D const& upperLeftCorner,
        Index2D const& lowerRightCorner)
    {
      Type2 result = static_cast<Type2>(0);
      for(int row = upperLeftCorner.getRow();
          row < lowerRightCorner.getRow();
          ++row) {
        typename Array2D<Type>::const_iterator rowBegin = array0.rowBegin(row);
        result += std::accumulate(
          rowBegin + upperLeftCorner.getColumn(),
          rowBegin + lowerRightCorner.getColumn(),
          static_cast<Type2>(0));
      }
      return result;
    }


    // This function returns an array made up of only those elements
    // of dataArray that correspond to indices in indexArray.
    template <class Type, class IntegralType>
    Array1D<Type>
    take(Array1D<Type> const& dataArray,
         Array1D<IntegralType> const& indexArray)
    {
      Array1D<Type> resultArray(indexArray.size());
      for(unsigned int ii = 0; ii < indexArray.size(); ++ii) {
        resultArray[ii] = dataArray[indexArray[ii]];
      }
      return resultArray;
    }


    // This function returns an array made up of only those elements
    // of dataArray that correspond to indices in indexArray.
    template <class Type, class IntegralType>
    Array1D<Type>
    take(Array2D<Type> const& dataArray,
         Array1D<IntegralType> const& indexArray)
    {
      return take(dataArray.ravel(), indexArray);
    }


    // This function works just like take(Array2D const&, Array1D
    // const&), with the exception that the input array is not
    // flattened, and entire rows (or columns) are selected.
    template <class Type, class IntegralType>
    Array2D<Type>
    take(Array2D<Type> const& dataArray,
         Array1D<IntegralType> const& indexArray,
         unsigned int axis)
    {
      Array2D<Type> resultArray;
      switch(axis) {
      case 0:
        resultArray.reinit(indexArray.size(), dataArray.columns());
        for(unsigned int ii = 0; ii < indexArray.size(); ++ii) {
          unsigned int jj = indexArray[ii];
          std::copy(dataArray.rowBegin(jj), dataArray.rowEnd(jj),
                    resultArray.rowBegin(ii));
        }
        break;
      default:
        resultArray.reinit(dataArray.rows(), indexArray.size());
        for(unsigned int ii = 0; ii < dataArray.rows(); ++ii) {
          typename Array2D<Type>::iterator resultIterator =
            resultArray.rowBegin(ii);
          typename Array2D<Type>::const_iterator dataRowBegin =
            dataArray.rowBegin(ii);
          for(unsigned int jj = 0; jj < indexArray.size(); ++jj) {
            *resultIterator = *(dataRowBegin + indexArray[jj]);
            ++resultIterator;
          }
        }
        break;
      }
      return resultArray;
    }


    // This function computes the standard deviation of a group of
    // scalar samples.
    template <class Type0, class Type1>
    inline Type1
    variance(Array1D<Type0> const& array0)
    {
      Type1 meanValue = mean<Type1>(array0);
      Type1 accumulator = 0;
      for(size_t index0 = 0; index0 < array0.size(); ++index0) {
        Type1 offset = static_cast<Type1>(array0[index0]) - meanValue;
        accumulator += offset * offset;
      }
      return accumulator / static_cast<Type1>(array0.size());
    }


    // This function returns an Array1D of the specified size and type
    // in which the value of every element is zero.
    template <class Type>
    inline Array1D<Type>
    zeros(size_t size)
    {
      Array1D<Type> result(size);
      result = 0;
      return result;
    }


    // This function returns an Array2D of the specified size and type
    // in which the value of every element is zero.
    template <class Type>
    inline Array2D<Type>
    zeros(size_t rows, size_t columns)
    {
      Array2D<Type> result(rows, columns);
      result = 0;
      return result;
    }


    // This function returns an Array3D of the specified size and type
    // in which the value of every element is zero.
    template <class Type>
    inline Array3D<Type>
    zeros(size_t shape0, size_t shape1, size_t shape2)
    {
      Array3D<Type> result(shape0, shape1, shape2);
      result = 0;
      return result;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_UTILITIES_IMPL_HH */
