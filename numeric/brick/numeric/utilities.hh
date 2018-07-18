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

#ifndef BRICK_NUMERIC_UTILITIES_HH
#define BRICK_NUMERIC_UTILITIES_HH

#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/array3D.hh>
#include <brick/numeric/index2D.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/numeric/vector3D.hh>

namespace brick {

  namespace numeric {

    /**
     * This function returns an array of the same size and element type
     * as its input argument, in which each element is set to the
     * absolute value of the corresponding element of the input array.
     *
     * @param array0 The elements of the returned array will take on the
     * absolute value of the elements of this array.
     *
     * @return The return value is an array in which each element
     * contains the absolute value of the corresponding element in the
     * input array.
     */
    template <class Type>
    Array1D<Type>
    abs(Array1D<Type> const& array0);


    /**
     * This function returns an array of the same size and element type
     * as its input argument, in which each element is set to the
     * absolute value of the corresponding element of the input array.
     *
     * @param array0 The elements of the returned array will take on the
     * absolute value of the elements of this array.
     *
     * @return The return value is an array in which each element
     * contains the absolute value of the corresponding element in the
     * input array.
     */
    template <class Type>
    Array2D<Type>
    abs(Array2D<Type> const& array0);


    /**
     * This function returns true if each element of its argument is
     * false, and returns false otherwise.
     *
     * @param array0 The elements of this array will be tested.
     *
     * @return The return value is true if all of the elements of array0
     * are false.
     */
    template <class Type>
    bool
    allFalse(Array1D<Type> const& array0);


    /**
     * This function returns true if each element of its argument is
     * true, and returns false otherwise.
     *
     * @param array0 The elements of this array will be tested.
     *
     * @return The return value is true if all of the elements of array0
     * are true.
     */
    template <class Type>
    bool
    allTrue(Array1D<Type> const& array0);


    /**
     * This function returns true if any element of its argument is
     * true, and returns false otherwise.
     *
     * @param array0 The elements of this array will be tested.
     *
     * @return The return value is true if any of the elements of array0
     * are false.
     */
    template <class Type>
    inline bool
    anyFalse(Array1D<Type> const& array0);


    /**
     * This function returns true if any element of its argument is
     * true, and returns false otherwise.
     *
     * @param array0 The elements of this array will be tested.
     *
     * @return The return value is true if any of the elements of array0
     * are true.
     */
    template <class Type>
    inline bool
    anyTrue(Array1D<Type> const& array0);


    /**
     * This function returns the index of the largest element of its input
     * array.  This function is equivalent to the quantity
     *   (std::max_element(array0.begin(), array0.end()) - array0.begin())
     *
     * @param array0 The elements of this array will be evaluated, and
     * the index of the largest will be returned.
     *
     * @return The index of the largest element of the array.
     */
    template <class Type>
    inline size_t
    argmax(Array1D<Type> const& array0);


    /**
     * This function returns the index of the largest element of its input
     * sequence.  This function is equivalent to the quantity
     *   (std::max_element(beginIter, endIter) - beginIter)
     *
     * @param beginIter This argument marks the start of the input
     * sequence.
     *
     * @param endIter This argument marks (one element past) the end
     * of the input sequence.
     *
     * @return The index of the largest element of the sequence.
     */
    template <class IterType>
    inline size_t
    argmax(IterType beginIter, IterType endIter);


    /**
     * This function returns the index of the largest element of its
     * input array, where largeness is defined by the second argument.
     * This function is equivalent to the quantity
     *
     *   (std::max_element(array0.begin(), array0.end(), comparator)
     *    - array0.begin())
     *
     * NOTE: Read the argument description for comparator carefully.  It
     * is consistent with the standard library convention, but many
     * people find it to be counterintuitive.
     *
     * @param array0 The elements of this array will be evaluated, and
     * the index of the largest will be returned.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return The index of the largest element of the array.
     */
    template <class Type, class Functor>
    inline size_t
    argmax(Array1D<Type> const& array0, Functor comparator);


    /**
     * This function returns the index of the largest element of its
     * input sequence, where largeness is defined by the second argument.
     * This function is equivalent to the quantity
     *
     *   (std::max_element(beginIter, endIter, comparator)
     *    - array0.begin())
     *
     * NOTE: Read the argument description for comparator carefully.  It
     * is consistent with the standard library convention, but many
     * people find it to be counterintuitive.
     *
     * @param beginIter This argument marks the start of the input
     * sequence.
     *
     * @param endIter This argument marks (one element past) the end
     * of the input sequence.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return The index of the largest element of the array.
     */
    template <class IterType, class Functor>
    inline size_t
    argmax(IterType beginIter, IterType endIter, Functor comparator);


    /**
     * This function returns an Index2D instance indicating which is
     * the largest element of its input array.
     *
     * @param array0 The elements of this array will be evaluated, and
     * the index of the largest will be returned.
     *
     * @return The index of the largest element of the array.
     */
    template <class Type>
    inline Index2D
    argmax2D(Array2D<Type> const& array0);


    /**
     * This function returns an Index2D instance indicating which is
     * the largest element of its input array, where largeness is
     * defined by the second argument.
     *
     * NOTE: Read the argument description for comparator carefully.
     * It is consistent with the standard library convention, but many
     * people find it to be counterintuitive.
     *
     * @param array0 The elements of this array will be evaluated, and
     * the index of the largest will be returned.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return The index of the largest element of the array.
     */
    template <class Type, class Functor>
    inline Index2D
    argmax2D(Array2D<Type> const& array0, Functor comparator);


    /**
     * This function returns the index of the smallest element of its input
     * array.  This function is equivalent to
     *   (std::min_element(array0.begin(), array0.end()) - array0.begin());
     *
     * @param array0 The elements of this array will be evaluated, and
     * the index of the smallest will be returned.
     *
     * @return The index of the smallest element of the array.
     */
    template <class Type>
    inline size_t
    argmin(Array1D<Type> const& array0);


    /**
     * This function returns the index of the smallest element of its
     * input array, where smallness is defined by the second argument.
     * This function is equivalent to the quantity
     *
     *   (std::min_element(array0.begin(), array0.end(), comparator)
     *    - array0.begin())
     *
     * @param array0 The elements of this array will be evaluated, and
     * the index of the smallest will be returned.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return The index of the smallest element of the array.
     */
    template <class Type, class Functor>
    inline size_t
    argmin(Array1D<Type> const& array0, Functor comparator);


    /**
     * This function returns an Index2D instance indicating which is
     * the smallest element of its input array.
     *
     * @param array0 The elements of this array will be evaluated, and
     * the index of the smallest will be returned.
     *
     * @return The index of the smallest element of the array.
     */
    template <class Type>
    inline Index2D
    argmin2D(Array2D<Type> const& array0);


    /**
     * This function returns an Index2D instance indicating which is
     * the smallest element of its input array, where smallness is
     * defined by the second argument.
     *
     * @param array0 The elements of this array will be evaluated, and
     * the index of the smallest will be returned.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return The index of the smallest element of the array.
     */
    template <class Type, class Functor>
    inline Index2D
    argmin2D(Array2D<Type> const& array0, Functor comparator);


    /**
     * This function returns an array of indices, result, so that the
     * sequence (array0[result[0]], array0[result[1]],
     * array0[result[2]], ...) is sorted from smallest to largest using
     * operator<().
     *
     * @param array0 The of this array will be compared to each other in
     * order to establish the correct sequence of indices.
     *
     * @return An array of indices as described above.
     *
     * Note(xxx): provide a version that uses iterators.
     */
    template <class Type>
    Array1D<size_t>
    argsort(Array1D<Type> const& array0);


//   /**
//    * This function returns an array of indices, result, so that the
//    * sequence (array0[result[0]], array0[result[1]],
//    * array0[result[2]], ...) is sorted from smallest to largest using
//    * the supplied comparison operator.
//    *
//    * @param array0 The of this array will be compared to each other in
//    * order to establish the correct sequence of indices.
//    *
//    * @param comparator This argument is a functor which supports the
//    * std::binary_function<Type, Type, bool> interface, and returns
//    * true if its first argument is smaller than its second argument.
//    *
//    * @return An array of indices as described above.
//    */
//   template <class Type, class Functor>
//   Array1D<size_t>
//   argsort(Array1D<Type> const& array0, Functor comparator);


    /**
     * This function is just an alias for the function axisMaximum().
     *
     * @param array0 See documentation for axisMaximum().
     *
     * @param axis See documentation for axisMaximum().
     *
     * @return See documentation for axisMaximum().
     */
    template <class Type>
    Array1D<Type>
    axisMax(Array2D<Type> const& array0,
            size_t axis) {return axisMaximum(array0, axis);}

    /**
     * This function returns an Array1D in which each element has the
     * value of the largest element in one row or column of the input
     * Array2D.  The maximum is taken along the axis specified by the
     * second argument.  For more information see the documentation of
     * axisMaximum(const Array2D&, size_t, Functor).
     *
     * @param array0 This is the input Array2D instance.
     *
     * @param axis This argument indicates along which axis the maximum
     * should be computed. Axis 0 is aligned with the columns of the
     * array, while axis 1 is aligned with the rows of the array.
     *
     * @return The return value is an Array1D instance containing the
     * maximum element from each row or column of the input array.
     */
    template <class Type>
    inline Array1D<Type>
    axisMaximum(Array2D<Type> const& array0, size_t axis);

    /**
     * This function returns an Array1D in which each element has the
     * value of the largest element in one row or column of the input
     * Array2D, where largeness is defined by the third argument.  The
     * maximum is taken along the axis specified by the second argument.
     * For example:
     *
     *   Array2D<int> testArray("[[1, 3, 0], [5, 0, 2]]");
     *   std::cout << axisMaximum(testArray, 0, std::less<int>()) << ", "
     *             << axisMaximum(testArray, 1, std::less<int>()) << ".";
     *
     * will print "Array1D([5, 3, 2]), Array1D([3, 5])."
     *
     * NOTE: Read the argument description for comparator carefully.  It
     * is consistent with the standard library convention, but many
     * people find it to be counterintuitive.
     *
     * @param array0 This is the input Array2D instance.
     *
     * @param axis This argument indicates along which axis the maximum
     * should be computed. Axis 0 is aligned with the columns of the
     * array, while axis 1 is aligned with the rows of the array.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return The return value is an Array1D instance containing the
     * maximum element from each row or column of the input array.
     */
    template <class Type, class Functor>
    Array1D<Type>
    axisMaximum(Array2D<Type> const& array0, size_t axis, Functor comparator);

    /**
     * This function is just an alias for the function axisMinimum().
     *
     * @param array0 See documentation for axisMinimum().
     *
     * @param axis See documentation for axisMinimum().
     *
     * @return See documentation for axisMinimum().
     */
    template <class Type>
    inline Array1D<Type>
    axisMin(Array2D<Type> const& array0,
            size_t axis) {return axisMinimum(array0, axis);}

    /**
     * This function returns an Array1D in which each element has the
     * value of the smallest element in one row or column of the input
     * Array2D.  The minimum is taken along the axis specified by the
     * second argument.  For more information see the documentation of
     * axisMinimum(const Array2D&, size_t, Functor).
     *
     * @param array0 This is the input Array2D instance.
     *
     * @param axis This argument indicates along which axis the minimum
     * should be computed. Axis 0 is aligned with the columns of the
     * array, while axis 1 is aligned with the rows of the array.
     *
     * @return The return value is an Array1D instance containing the
     * minimum element from each row or column of the input array.
     */
    template <class Type>
    inline Array1D<Type>
    axisMinimum(Array2D<Type> const& array0, size_t axis);

    /**
     * This function returns an Array1D in which each element has the
     * value of the smallest element in one row or column of the input
     * Array2D, where smallness is defined by the return value of the
     * third argument.  The minimum is taken along the axis specified
     * by the second argument.  For example:
     *
     *   Array2D<int> testArray("[[1, 3, 0], [5, 0, 2]]");
     *   std::cout << axisMinimum(testArray, 0, std::less<int>()) << ", "
     *             << axisMinimum(testArray, 1, std::less<int>()) << ".";
     *
     * will print "Array1D([1, 0, 0]), Array1D([0, 0])."
     *
     * NOTE: Read the argument description for comparator carefully.  It
     * is consistent with the standard library convention, but many
     * people find it to be counterintuitive.
     *
     * @param array0 This is the input Array2D instance.
     *
     * @param axis This argument indicates along which axis the minimum
     * should be computed. Axis 0 is aligned with the columns of the
     * array, while axis 1 is aligned with the rows of the array.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return The return value is an Array1D instance containing the
     * minimum element from each row or column of the input array.
     */
    template <class Type, class Functor>
    Array1D<Type>
    axisMinimum(Array2D<Type> const& array0, size_t axis, Functor comparator);


    /**
     * This function returns an Array1D in which each element has the
     * sum of one row or column of the input Array2D.  The sum is taken
     * along the axis specified by the second argument.  For Example:
     *
     *   Array2D<int> testArray("[[1, 3, 0], [5, 0, 2]]");
     *   std::cout << axisSum<int>(testArray, 0) << ", "
     *             << axisSum<int>(testArray, 1) << ".";
     *
     * will print "Array1D([6, 3, 2]), Array1D([4, 7])."
     *
     * The summation is done using values of type
     * NumericTraits<Type>::SumType, where Type is the type of the
     * elements in the input array.
     *
     * AxisSum permits the user to control the precision of the
     * summation by specifying the type which will be used to do the
     * computation using the second template argument.
     *
     * @param array0 This is the input Array2D instance.
     *
     * @param axis This argument indicates along which axis the sum
     * should be computed. Axis 0 is aligned with the columns of the
     * array, while axis 1 is aligned with the rows of the array.
     *
     * @return The return value is an Array1D instance containing the
     * sums of the rows or columns of the input array.
     */
    template <class ResultType, class Type>
    inline Array1D<ResultType>
    axisSum(Array2D<Type> const& array0, size_t axis);


    /**
     * This function returns an Array1D in which each element has the
     * sum of one row or column of the input Array2D.  The sum is taken
     * along the axis specified by the second argument, as described for
     * axisSum(const Array2D&, size_t), above, and addition is performed
     * using a user supplied binary operator.  For each row or column,
     * the sum is initialized to the value specified by initialValue
     * before the summation is performed.  After this initialization,
     * each row- or column- sum is updated according to the equation:
     *
     *   sum = adder(sum, nextElement)
     *
     * for each element of the row or column.  For more information, see
     * the documentation of other axisSum() versions.
     *
     * AxisSum also permits the user to control the precision of the
     * summation by specifying (via the second template argument) the
     * type that will be used to do the computation.
     *
     * @param array0 This is the input Array2D instance.
     *
     * @param axis This argument indicates along which axis the sum
     * should be computed. Axis 0 is aligned with the columns of the
     * array, while axis 1 is aligned with the rows of the array.
     *
     * @param initialValue This argument specifies the value to which
     * each row/column sum will be initialized prior to accumulating the
     * elements of the array.
     *
     * @param adder This argument is a functor conforming to the
     * std::binary_function interface which will be used to add element
     * values.
     *
     * @return The return value is an Array1D instance containing the
     * sums of the rows or columns of the input array.
     */
    template <class ResultType, class Type, class Functor>
    Array1D<ResultType>
    axisSum(Array2D<Type> const& array0, size_t axis,
            ResultType const& initialValue, Functor adder);


    /**
     * This function swaps values between its two arguments so that
     * the value of each element of corner0 is less than or equal to
     * the corresponding value of corner1.  After calling this
     * function, you can safely write code that depends on the
     * relative positioning of the two corners.
     *
     * @param corner0 This argument is the corner that should have
     * smaller values.
     *
     * @param corner1 This argument is the corner that should have
     * larger values.
     *
     * @return The return value is true if -- after the call to
     * cleanupCorners() -- the area of the rectangle defined by
     * corner0 and corner1 is non-zero.  That is, the return value is
     * false if and only if ((corner0.getX() == corner1.getX()) ||
     * (corner0.getY() == corner1.getY())).
     */
    template <class Type>
    bool
    cleanupCorners(Vector2D<Type>& corner0, Vector2D<Type>& corner1);


    /**
     * columnIndices(rows, columns): Returns an Array2D in which each
     * element contains the index of its column.  The element type of the
     * returned array depends on the type_tag argument.  For example:
     *
     *   columnIndices<int>(2, 3);
     *
     * will return:
     *
     *   Array2D<int>([[0, 1, 2],
     *                 [0, 1, 2]]);
     *
     *
     * @param rows This argument specifies the number of rows in the
     * returned array.
     *
     * @param columns This argument specifies the number of columns in
     * the returned array.
     *
     * @return The return value is an Array2D of the specified size and
     * type in which all elements in the first column have value
     * initialized to 0, all elements in the second column have value
     * initialized to 1, and so on.
     */
    template <class Type>
    Array2D<Type>
    columnIndices(size_t rows, size_t columns);


    /**
     * This function selects those elements of an input Array1D which
     * correspond to "true" values of a mask Array1D, and returns an
     * Array1D containing only those elements.  For example:
     *
     *   Array1D<int> a = compress(
     *     Array1D<bool>("1, 1, 0, 1, 0]"),
     *     Array1D<int>("1, 2, 3, 4, 5"));
     *   std::cout << a;
     *
     * will print:
     *
     *   "Array1D([1, 2, 4])"
     *
     * It is an error for the two array arguments to have different sizes.
     *
     * @param condition This argument is the mask which selects which
     * elements will be included in the output array.  Elements of
     * argument this array which evaluate to true will cause the
     * corresponding elements of the input array to be included in the
     * output.
     *
     * @param input Elements of this argument will be included in the
     * output array if the corresponding elements of the condition array
     * evaluate to true.  This array must have the same number of
     * elements as the condition array.
     *
     * @return The return value is an Array1D instance containing only
     * those elements of the input array which correspond to true
     * elements of the condition array.
     */
    template <class Type0, class Type1>
    inline Array1D<Type1>
    compress(Array1D<Type0> const& condition, Array1D<Type1> const& input);

    /**
     * This function behaves in exactly the same way as compress(const
     * Array1D&, const Array1D&), above, but it permits the user to
     * specify the number of true elements in the condition array.  This
     * eliminates one pass through the condition array, slightly
     * speeding up the computation.  Use this version with caution,
     * since specifying an incorrect count can crash your program.
     *
     * @param condition This argument is the mask which selects which
     * elements will be included in the output array.  Elements of
     * argument this array which evaluate to true will cause the
     * corresponding elements of the input array to be included in the
     * output.
     *
     * @param input Elements of this argument will be included in the
     * output array if the corresponding elements of the condition array
     * evaluate to true.  This array must have the same number of
     * elements as the condition array.
     *
     * @param numTrue This argument specifies the number of "true"
     * elements in the condition array.
     *
     * @return The return value is an Array1D instance containing only
     * those elements of the input array which correspond to true
     * elements of the condition array.
     */
    template <class Type0, class Type1>
    Array1D<Type1>
    compress(Array1D<Type0> const& condition, Array1D<Type1> const& input,
             size_t numTrue);


    /**
     * This element counts the number of elements of the input array
     * which evaluate to true, and returns that number.
     *
     * @param array0 This argument is the array of elements to be
     * counted.
     *
     * @return The return value is the total number of elements in x
     * which, when cast to bool, evaluate to true.
     */
    template <class Type>
    inline size_t
    count(Array1D<Type> const& array0);

    /**
     * This function computes the cross product of two Vector3D
     * instances.  That is, it computes a Vector3D which lies
     * perpendicular to each of the arguments and has magnitude equal
     * to the product of the magnitudes of the arguments times the sine
     * of the angle between them.
     *
     * @param vector0 This is the first argument of the cross product.
     *
     * @param vector1 This is the first argument of the cross product.
     *
     * @return The cross product of the two arguments.
     */
    template <class Type>
    inline Vector3D<Type>
    cross(Vector3D<Type> const& vector0, Vector3D<Type> const& vector1);

    /**
     * This function computes the inner product of two input arrays.
     * The computation is done using the ProductType specified by the
     * third template argument.
     *
     * @param array0 A 1D array, the first argument of the dot product.
     *
     * @param array1 A 1D array, the second argument of the dot product.
     *
     * @return The inner product of the two arguments.
     */
    template <class Type2, class Type1, class Type0>
    inline Type2
    dot(Array1D<Type0> const& array0, Array1D<Type1> const& array1);


    /**
     * This function computes the inner product of two Vector2D instances.
     *
     * @param vector0 A Vector2D, the first argument of the dot product.
     *
     * @param vector1 A Vector2D, the second argument of the dot product.
     *
     * @return The inner product of the two Vector2D arguments.
     */
    template <class Type1, class Type0>
    inline Type1
    dot(Vector2D<Type0> const& vector0, Vector2D<Type0> const& vector1);


    /**
     * This function computes the inner product of two Vector3D instances.
     *
     * @param vector0 A Vector3D, the first argument of the dot product.
     *
     * @param vector1 A Vector3D, the second argument of the dot product.
     *
     * @return The inner product of the two Vector3D arguments.
     */
    template <class Type1, class Type0>
    inline Type1
    dot(Vector3D<Type0> const& vector0, Vector3D<Type0> const& vector1);


    /**
     * This function computes from a vector, x, the matrix, X, such that
     * the matrix product A * x is equal to the matrix product X *
     * vec(A).  That is, if
     *
     *   X = equivalentMatrix(x, n);
     *
     * then
     *
     *   A * x == X * vec(A)
     *
     * where x is a vector of m elements, and A is a matrix with n rows
     * and m columns.
     *
     * @param vector0 This argument represents the vector x.
     *
     * @param rowsInMatrix This argument describes the matrix A.
     *
     * @return The return value is an Array2D instance describing the
     * matrix X.
     */
    template <class Type>
    Array2D<Type>
    equivalentMatrix(Array1D<Type> const& vector0, size_t rowsInMatrix);


    /**
     * This function returns the centroid of a 1D array.  That is, it
     * returns the value
     *
     * @code
     *   centroid = sum_i(f[i] * i) / sum_i(f[i])
     * @endcode
     *
     * where f[i] is the i^th element of the input sequence.
     *
     * @param signal This argument is the input sequence.
     *
     * @return The return value is the computed centroid.
     */
    template <class FloatType, class Type>
    FloatType
    getCentroid(Array1D<Type> const& signal);


    /**
     * This function returns the centroid of a 1D array.  That is, it
     * returns the value
     *
     * @code
     *   centroid = sum_i(f[i] * i) / sum_i(f[i])
     * @endcode
     *
     * where f[i] is the i^th element of the input sequence.
     *
     * @param signal This argument is the input sequence.
     *
     * @return The return value is the computed centroid.
     */
    template <class FloatType, class IterType>
    FloatType
    getCentroid(IterType beginIter, IterType endIter);


    /**
     * This function estimates the mean and covariance of a set of
     * vectors, which are represented by the rows (or columns) of the
     * input 2D array.  Estimated mean and covariance are returned via
     * reference arguments.
     *
     * @param sampleArray This argument is the array of sample vectors.
     * If argument majorAxis is set to 0, then each row of sampleArray
     * represents a sample.  If argument majorAxis is nonzero, then each
     * column of sampleArray represents a sample.
     *
     * @param mean The estimated mean vector is returned by reference
     * using this argument.
     *
     * @param covariance The estimated covariance matrix is returned by
     * reference using this argument.
     *
     * @param majorAxis This argument is normally set to 0 or 1, and
     * specifies the arrangement of samples in sampleArray as described
     * above.  If majorAxis is set to a value which is not equal to 0 or
     * 1, computation will proceed as if it had been set to 1.
     */
    template <class Type>
    void
    getMeanAndCovariance(Array2D<Type> const& sampleArray,
                         Array1D<double>& mean,
                         Array2D<double>& covariance,
                         size_t majorAxis=0);


    /**
     * This function estimates the mean and variance of a sequence of
     * scalar values.
     *
     * Template argument Iter is the type of iterator used to pass in
     * the sequence of scalars.
     *
     * Template argument Type defines the numeric type that will be
     * used during computation of the mean and variance.
     *
     * @param beginIter This iterator, and the next, define sequence
     * of scalars.
     *
     * @param endIter This iterator, and the previous, define sequence
     * of scalars.
     *
     * @param mean The estimated mean vector is returned by reference
     * using this argument.
     *
     * @param variance The estimated variance is returned by reference
     * using this argument.
     */
    template <class Iter, class Type>
    void
    getMeanAndVariance(Iter beginIter, Iter endIter,
                       Type& mean, Type& variance);


#if 0 /* Not implemented yet. */
    /**
     * This function is not currently implemented.  If it were, it
     * would estimate the mean and variance of a sequence of scalar
     * values, discarding values that appear to be outliers.
     *
     * Template argument Iter is the type of iterator used to pass in
     * the sequence of scalars.
     *
     * Template argument Type defines the numeric type that will be
     * used during computation of the mean and variance.
     *
     * @param beginIter This iterator, and the next, define sequence
     * of scalars.
     *
     * @param endIter This iterator, and the previous, define sequence
     * of scalars.
     *
     * @param inlierProportion Ths argument indicates what proportion
     * of the input values are expected to be inliers.  Setting this
     * to 1.0 (or greater) means that all input values are expected to
     * be valid.  Setting this to 0.0 (or a negative number) indicates
     * that no input values are expected to be valid, and is an error.
     *
     * @param mean The estimated mean vector is returned by reference
     * using this argument.
     *
     * @param variance The estimated variance is returned by reference
     * using this argument.
     */
    template <class Iter, class Type>
    void
    getMeanAndVarianceRobust(Iter beginIter, Iter endIter,
                             Type& mean, Type& variance);
#endif /* #if 0 */


    /**
     * This function returns an Array2D instance with the specified
     * shape in which the elements on the diagonal are set to 1, and all
     * other elements are set to 0.
     *
     * @param rows This argument specifies the number of rows in the
     * returned Array2D instance.
     *
     * @param columns This argument specifies the number of columns in
     * the returned Array2D instance.
     *
     * @return An Array2D in each element on the diagonal is set to 1,
     * and all other elements are set to 0.
     */
    template <class Type>
    inline Array2D<Type>
    identity(size_t rows, size_t columns);


    /**
     * This function returns an array in which each element is the
     * natural logarithm of the corresponding element of its input.
     *
     * @param array0 This argument is an Array1D with arguments in a
     * common floating point type such as double or float.
     *
     * @return The return value is an array of natural log values.
     */
    template <class Type>
    Array1D<Type>
    ln(Array1D<Type> const& array0);


    /**
     * This function returns an array in which each element is the
     * natural logarithm of the corresponding element of its input.
     *
     * @param array0 This argument is an Array2D with arguments in a
     * common floating point type such as double or float.
     *
     * @return The return value is an array of natural log values.
     */
    template <class Type>
    Array2D<Type>
    ln(Array2D<Type> const& array0);


    /**
     * This function returns an Array2D instance in which the value of
     * each element is the logical not of the corresponding element of
     * the input array.  For example:
     *
     *   Array2D<int> a = logical_not(Array2D<int>("[[1, 2], [0, -5]]");
     *   std::cout << a;
     *
     * will print:
     *
     *   "Array2D([[0, 0], [1, 0]])"
     *
     * @param array0 This argument is an Array2D with elements of a type
     * which can be cast to bool.
     *
     * @return The return value is an Array2D instance of the same shape
     * and element type as the input array, in which each element is set
     * to the logical not of the corresponding element of the input
     * array.
     */
    template <class Type>
    Array2D<Type>
    logicalNot(Array2D<Type> const& array0);

    /**
     * This function computes the magnitude of its input argument.  That
     * is, it computes the square root of the sum of the squares of the
     * elements of its input argument.
     *
     * @param array0 An Array1D whose magnitude is to be calculated.
     *
     * @return The magnitude of the input argument.
     */
    template <class Type1, class Type0>
    inline Type1
    magnitude(Array1D<Type0> const& array0);

    /**
     * This function computes the magnitude of its input argument.  That
     * is, it computes the square root of the sum of the squares of the
     * elements of its input argument.
     *
     * @param vector0 A Vector2D whose magnitude is to be calculated.
     *
     * @return The magnitude of the input argument.
     */
    template <class Type1, class Type0>
    inline Type1
    magnitude(Vector2D<Type0> const& vector0);

    /**
     * This function computes the magnitude of its input argument.  That
     * is, it computes the square root of the sum of the squares of the
     * elements of its input argument.
     *
     * @param vector0 A Vector3D whose magnitude is to be calculated.
     *
     * @return The magnitude of the input argument.
     */
    template <class Type1, class Type0>
    inline Type1
    magnitude(Vector3D<Type0> const& vector0);


    /**
     * This function computes the square of the magnitude of its input
     * argument.  That is, it computes the sum of the squares of the
     * elements of its input argument.
     *
     * @param array0 An Array1D whose squared magnitude is to be
     * calculated.
     *
     * @return The squared magnitude of the input argument.
     */
    template <class Type1, class Type0>
    inline Type1
    magnitudeSquared(Array1D<Type0> const& array0);


    /**
     * This function computes the square of the magnitude of its input
     * argument.  That is, it computes the sum of the squares of the
     * elements of its input argument.
     *
     * @param vector0 A Vector2D whose squared magnitude is to be
     * calculated.
     *
     * @return The squared magnitude of the input argument.
     */
    template <class Type1, class Type0>
    inline Type1
    magnitudeSquared(Vector2D<Type0> const& vector0);


    /**
     * This function computes the square of the magnitude of its input
     * argument.  That is, it computes the sum of the squares of the
     * elements of its input argument.
     *
     * @param vector0 A Vector3D whose squared magnitude is to be
     * calculated.
     *
     * @return The squared magnitude of the input argument.
     */
    template <class Type1, class Type0>
    inline Type1
    magnitudeSquared(Vector3D<Type0> const& vector0);


    /**
     * This function computes a vector * matrix product.  Assuming the
     * first argument, vector0, represents a column vector and the
     * second argument, matrix0, represents a 2D array, then this function
     * computes the product:
     *
     *  vector0' * matrix0
     *
     * Where the single quote indicates transpose.
     *
     * The element type of the return value is set explicitly using
     * the third template argument.
     *
     * @param vector0 The first operand for the multiplication.
     *
     * @param matrix0 The second operand for the multiplication.
     *
     * @return A 1D array which is the matrix product of vector0 and
     * matrix0.  The element type of this array is determined by
     * NumericTraits<Type>::ProductType for the element type of the
     * arguments.
     */
    template <class Type2, class Type1, class Type0>
    Array1D<Type2>
    matrixMultiply(Array1D<Type0> const& vector0,
                   Array2D<Type1> const& matrix0);


    /**
     * This function computes a matrix * vector product.  Assuming the
     * first argument, matrix0, represents a matrix, and the and the
     * second argument, vector0, represents a column vector, then this
     * function computes the product:
     *
     *  matrix0 * vector0
     *
     * The element type of the return value is set explicitly using a
     * third template argument.
     *
     * @param matrix0 The first operand for the multiplication.
     *
     * @param vector0 The second operand for the multiplication.
     *
     * @return A 1D array which is the matrix product of matrix0 and
     * vector0.  The element type of this array is determined by
     * NumericTraits<Type>::ProductType for the element type of the
     * arguments.
     */
    template <class Type2, class Type1, class Type0>
    Array1D<Type2>
    matrixMultiply(Array2D<Type0> const& matrix0,
                   Array1D<Type1> const& vector0);


    /**
     * This function computes a matrix * matrix product.  That is,
     * the elements of the resulting array are the dot products of the
     * rows of the first argument with the columns of the second argument.
     *
     * The element type of the return value is set explicitly using the
     * third template argument.
     *
     * @param matrix0 The first operand for the multiplication.
     *
     * @param matrix1 The second operand for the multiplication.
     *
     * @return A 2D array which is the matrix product of matrix0 and
     * matrix1.  The element type of this array is determined by
     * NumericTraits<Type>::ProductType for the element type of the
     * arguments.
     */
    template <class Type2, class Type1, class Type0>
    Array2D<Type2>
    matrixMultiply(Array2D<Type0> const& matrix0,
                   Array2D<Type1> const& matrix1);


    // The next function is commented out because some environments use
    // #define directives for min() and max(), which don't obey
    // namespaces, and makes the code not compile if we define our own
    // min/max.
    //
    //
    // /**
    // * This function is an alias for the function maximum(const
    // * Array1D&) below.
    // *
    // * @param array0 See documentation for maximum(Array1D const&).
    // *
    // * @return See documentation for maximum(Array1D const&).
    // */
    // template <class Type>
    // inline Type
    // max(Array1D<Type> const& array0) {return maximum(array0);}


    /**
     * This function returns a copy of the largest element in the input
     * Array1D instance.
     *
     * @param array0 This argument is the input array.
     *
     * @return A copy of an element of the input array such that no
     * other element of the input array is larger.
     */
    template <class Type>
    inline Type
    maximum(Array1D<Type> const& array0);


    /**
     * This function returns a copy of the largest element in the input
     * Array1D instance, where largeness is defined by the return value
     * of the second argument.
     *
     * NOTE: Read the argument description for comparator carefully.  It
     * is consistent with the standard library convention, but many
     * people find it to be counterintuitive.
     *
     * @param array0 This argument is the input array.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return A copy of an element of the input array such that no
     * other element of the input array is larger, where largeness is
     * defined by comparator.
     */
    template <class Type, class Functor>
    Type
    maximum(Array1D<Type> const& array0, Functor comparator);


    /**
     * This function returns a copy of the largest element in the input
     * Array2D instance.
     *
     * @param array0 This argument is the input array.
     *
     * @return A copy of an element of the input array such that no
     * other element of the input array is larger.
     */
    template <class Type>
    inline Type
    maximum(Array2D<Type> const& array0);


    /**
     * This function returns a copy of the largest element in the input
     * Array2D instance, where largeness is defined by the return value
     * of the second argument.
     *
     * NOTE: Read the argument description for comparator carefully.  It
     * is consistent with the standard library convention, but many
     * people find it to be counterintuitive.
     *
     * @param array0 This argument is the input array.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return A copy of an element of the input array such that no
     * other element of the input array is larger, where largeness is
     * defined by comparator.
     */
    template <class Type, class Functor>
    Type
    maximum(Array2D<Type> const& array0, Functor comparator);


    /**
     * This function computes the average value, or geometric mean, of
     * the elements of its argument.  The computation is performed using
     * the precision specified by NumericTraits<Type>::SumType, but is
     * returned as an instance of Type.
     *
     * @param array0 This argument is the input array, from the elements
     * of which the mean will be computed.
     *
     * @return The return value is the mean of the elements of array0.
     */
    template <class Type, class Iterator>
    Type
    mean(Iterator beginIter, Iterator endIter);


    /**
     * This function computes the average value, or geometric mean, of
     * the elements of its argument.  The computation is performed
     * using the type specified by the second template argument.
     *
     * @param array0 This argument is the input array, from the elements
     * of which the mean will be computed.
     *
     * @return The return value is the mean of the elements of array0.
     */
    template <class Type1, class Type0>
    inline Type1
    mean(Array1D<Type0> const& array0);


    // The next function is commented out because some environments use
    // #define directives for min() and max(), which don't obey
    // namespaces, and makes the code not compile if we define our own
    // min/max.
    //
    //
    // /**
    // * This function is an alias for the function minimum(const
    // * Array1D&) below.
    // *
    // * @param array0 See documentation for minimum(Array1D const&).
    // *
    // * @return See documentation for minimum(Array1D const&).
    // */
    // template <class Type>
    // inline Type
    // min(Array1D<Type> const& array0) {return minimum(array0);}


    /**
     * This function returns a copy of the smallest element in the input
     * Array1D instance.
     *
     * @param array0 This argument is the input array.
     *
     * @return A copy of an element of the input array such that no
     * other element of the input array is smaller.
     */
    template <class Type>
    inline Type
    minimum(Array1D<Type> const& array0);


    /**
     * This function returns a copy of the smallest element in the input
     * Array1D instance, where largeness is defined by the return value
     * of the second argument.
     *
     * NOTE: Read the argument description for comparator carefully.  It
     * is consistent with the standard library convention, but many
     * people find it to be counterintuitive.
     *
     * @param array0 This argument is the input array.
     *
     * @param comparator This argument is a functor which supports the
     * std::binary_function<Type, Type, bool> interface, and returns
     * true if its first argument is smaller than its second argument.
     *
     * @return A copy of an element of the input array such that no
     * other element of the input array is smaller, where largeness is
     * defined by comparator.
     */
    template <class Type, class Functor>
    Type
    minimum(Array1D<Type> const& array0, Functor comparator);


    /**
     * This function uses the iterative method of Newton and Raphson to
     * search for a zero crossing of the supplied functor.  If convergence
     * fails, it throws a ValueException (blaming the argument for not
     * having a sufficiently findable zero crossing).
     *
     * @param startPoint This argument specifies the point at which to
     * start the search.
     *
     * @param objectiveFunction This argument is an instance of the
     * functor type.  The functor must accept a single argument of type
     * double, and return a double, and must provide a method
     * derivative(double), which returns the first derivative.
     *
     * @param epsilon This argument specifies the convergence tolerance.
     * Iteration will stop when the algorithm finds a value which, when
     * passed to objectiveFunction, results in a return value with
     * magnitude less than epsilon.
     *
     * @param maxIterations This argument specifies how many iterations
     * should be permitted before giving up the search (& throwing an
     * exception).
     *
     * @return The return value is a value which, when passed to
     * objectiveFunction, results in a return value with magnitude less
     * than epsilon.
     */
    template <class FUNCTOR>
    double newtonRaphson(double startPoint, FUNCTOR objectiveFunction,
                         double epsilon, size_t maxIterations);


    /**
     * This function computes the normalized correlation of two Array1D
     * arguments.  This is equivalent to (but slightly more efficient
     * than) independently normalizing the signal in each Array1D to
     * zero mean and unit variance, and then computing the dot product
     * of the two normalized Array1D elements.  The computation is
     * carried out in double precision floating point.  For more
     * information on normalized correlation, which is often called
     * "correlation coefficient," please refer to R. O. Duda and
     * P. E. Hart, Pattern Classification and Scene Analysis, Wiley, NY,
     * 1973.
     *
     * The computation is carried out using the type specified by the
     * second template argument.
     *
     * @param signal0 This argument is an Array1D instance containing the
     * first of the signals to be correlated.
     *
     * @param signal1 This argument is an Array1D instance containing the
     * second of the signals to be correlated.
     *
     * @return The return value is the normalized correlation of the two
     * signals.
     */
    template <class Type2, class Type>
    Type2
    normalizedCorrelation(Array1D<Type> const& signal0,
                          Array1D<Type> const& signal1);


    /**
     * This function returns an Array1D of the specified size and type
     * in which the value of every element is initialized to 1.
     *
     * @param size This argument specifies the number of elements in the
     * returned array.
     *
     * @return An Array1D instance in which each element is set to zero.
     */
    template <class Type>
    inline Array1D<Type>
    ones(int size);


    /**
     * This function returns an Array2D of the specified size and type
     * in which the value of every element is initialized to one.
     *
     * @param rows This argument specifies the number of rows in the
     * returned array.
     *
     * @param columns This argument specifies the number of columns in the
     * returned array.
     *
     * @return An Array2D instance in which each element is set to zero.
     */
    template <class Type>
    inline Array2D<Type>
    ones(int rows, int columns);


    /**
     * This function computes the outer product of two input Array1D
     * instances and allows the user to control which type is used to do
     * the calculation.  The computation is done using the ProductType
     * specified by the third template argument.  The result is
     * Array2D instance with values such that the element in the i-th
     * column and j-th row has a value equal to the product of the i-th
     * element of the second input argument and the j-th element of the
     * first input argument.
     *
     * @param array0 A 1D array, the first argument of the outer product.
     *
     * @param array1 A 1D array, the second argument of the outer product.
     *
     * @return The outer product of the two arguments.
     */
    template <class Type2, class Type0, class Type1>
    Array2D<Type2>
    outerProduct(Array1D<Type0> const& array0, Array1D<Type1> const& array1);


    /**
     * This function returns an Array1D in which the first element has
     * value equal to argument "start," and each subsequent element has
     * value equal to the previous element plus argument "stride."  The
     * length of the array is equal to largest integer less than the
     * quantity ((stop - start) / stride).  It is an error for stride to
     * have a different sign than (stop - start).  The element type of
     * the returned array is determined by the type of the arguments.
     * For example,
     *
     *  range(2, 7, 1);
     *
     * will return
     *
     *  Array1D([2, 3, 4, 5, 6]);
     *
     * and
     *
     *  range(2.0, -6.0, -1.5);
     *
     * will return
     *
     *  Array1D([2.0, 0.5, -1.0, -2.5, -4.0, -5.5]);
     *
     * @param start This argument specifies the value of the first
     * element of the returned Array1D instance.
     *
     * @param stop This argument specifies an exclusive bound for the
     * element values of the returned array.  If stride is less than
     * 0.0, then this will be an exclusive lower bound, else it will be
     * an exclusive upper bound.
     *
     * @param stride This argument specifies the increment which will be
     * added to each subsequent element in the output array.  It is an
     * error for stride to be equal to 0.
     *
     * @return The return value is an Array1D instance having regularly
     * spaced elements progressing from start toward stop.
     */
    template <class Type>
    Array1D<Type>
    range(Type start, Type stop, Type stride=1);


    /**
     * This function takes an Array2D argument and returns an Array1D
     * instance which references the same data.  That is, changing an
     * element of the Array1D instance will change the corresponding
     * element of the input Array2D instance.  The elements of the
     * returned Array1D are in native storage order for the input array.
     * Input arrays of type Array2D happen to have their elements in
     * row-major order.  Will this always be true?  Hmm...
     *
     * @param inputArray This argument is the Array2D instance from
     * which to take the data.
     *
     * @return An Array1D instance referencing the same data as the
     * input Array2D instance.
     */
    template <class Type>
    inline Array1D<Type>
    ravel(Array2D<Type>& inputArray) {return inputArray.ravel();}


    /**
     * This function takes a Array2D argument and returns a const
     * Array1D instance which references the same data.  It is
     * equivalent to ravel(Array2D const&), but for arrays.
     *
     * @param inputArray This argument is the Array2D instance from
     * which to take the data.
     *
     * @return An Array1D instance referencing the same data as the
     * input Array2D instance.
     */
    template <class Type>
    inline const Array1D<Type>
    ravel(Array2D<Type> const& inputArray) {return inputArray.ravel();}


    /**
     * This function computes the RMS (Root Mean Square) value of the
     * elements of its argument, and allows the user to specify the
     * precision with which the computation is carried out.  The
     * computation is performed using the type specified by the second
     * argument to the rms() function call.
     *
     * @param array0 This argument is the input array, from which the
     * RMS value will be computed.
     *
     * @param typeTag This argument specifies the type which
     * will be used to perform the computation, and also specifies the
     * return type of the function.
     *
     * @return The return is the RMS value of the elements of array0.
     */
    template <class Type1, class Type0>
    Type1
    rms(Array1D<Type0> const& array0);


    /**
     * rowIndices(rows, columns): Returns an Array2D in which each
     * element contains the index of its row.  For example:
     *
     *   rowIndices<int>(2, 3);
     *
     * will return:
     *
     *   Array2D<int>([[0, 0, 0],
     *                 [1, 1, 1]]);
     *
     * @param rows This argument specifies the number of rows in the
     * returned array.
     *
     * @param columns This argument specifies the number of columns in
     * the returned array.
     *
     * @return The return value is an Array2D of the specified size and
     * type in which all elements in the first row have value
     * initialized to 0, all elements in the second row have value
     * initialized to 1, and so on.
     */
    template <class Type>
    Array2D<Type>
    rowIndices(size_t rows, size_t columns);


    /**
     * This function returns true if the two arrays have the same shape,
     * false otherwise.  Two arrays are considered to have the same
     * shape if their sizes are equal along each axis.
     *
     * @param array0 This argument is the first of the two arrays to be
     * compared.
     *
     * @param array1 This argument is the second of the two arrays to be
     * compared.
     *
     * @return The return value is true if the two arrays have matching
     * shape, false otherwise.
     */
    template <class Type0, class Type1>
    inline bool
    shapeMatch(Array1D<Type0> const& array0, Array1D<Type1> const& array1);


    /**
     * This function returns true if the two arrays have the same shape,
     * false otherwise.  Two arrays are considered to have the same
     * shape if their sizes are equal along each axis.
     *
     * @param array0 This argument is the first of the two arrays to be
     * compared.
     *
     * @param array1 This argument is the second of the two arrays to be
     * compared.
     *
     * @return The return value is true if the two arrays have matching
     * shape, false otherwise.
     */
    template <class Type0, class Type1>
    inline bool
    shapeMatch(Array2D<Type0> const& array0, Array2D<Type1> const& array1);


    /**
     * This function returns true if the two arrays have the same shape,
     * false otherwise.  Two arrays are considered to have the same
     * shape if their sizes are equal along each axis.
     *
     * @param array0 This argument is the first of the two arrays to be
     * compared.
     *
     * @param array1 This argument is the second of the two arrays to be
     * compared.
     *
     * @return The return value is true if the two arrays have matching
     * shape, false otherwise.
     */
    template <class Type0, class Type1>
    inline bool
    shapeMatch(Array3D<Type0> const& array0, Array3D<Type1> const& array1);


    /**
     * skewSymmetric(x): Returns a skew symmetric matrix X such
     * that matrixMultiply(X, y) = cross(x, y)
     *
     * @param vector0 An Array1D instance, which must have exactly 3 elements.
     *
     * @return An Array2D instance as described above.
     */
    template <class Type>
    inline Array2D<Type>
    skewSymmetric(Array1D<Type> const& vector0);

#if 0
    /**
     * This function computes the real roots of the quadratic polynomial
     * c0*x^2 + c1*x + c2 = 0.
     *
     * If c0, c1, and c2 are arrays, then the returned roots will be
     * arrays of the same dimension, and will contain the roots of an
     * array of scalar quadratic equations.
     *
     * Note that the two well known versions of the quadratic formula:
     *
     *   x = (-c1 +|- sqrt(c1**2 - 4c0*c2)) / (2*c0)
     *
     * and
     *
     *   x = 2*c2 / (-c1 +|- sqrt(c1**2 - 4*c0*c2))
     *
     * both tend to be inaccurate when c0 and/or c2 are small, since
     * then the quantity (-c1 +|- sqrt(c1**2 - 4*c0*c2)) gets very
     * small and loses significance.  Instead we use the form
     * advocated by Press and Flannery (Numerical Recipes):
     *
     *   q = -(1/2)(c1 + sgn(c1)*sqrt(c1**2 - 4*c1*c2))
     *   x1 = q/c0, x2 = c2/q
     *
     * This is just the same as using both of the first two forms, each
     * one for the root at which it is most numerically stable.
     *
     * If checkValidity is set to zero, valid will always be set to
     * (scalar) one, and no attempt will be made to detect whether
     * quadratics have real roots.  This will make your code run
     * faster, but could cause problems if your quadratics don't
     * always have real roots.
     *
     * @param c0 This argument is the quadratic coefficient of the
     * polynomial.
     *
     * @param c1 This argument is the linear coefficient of the
     * polynomial.
     *
     * @param c2 This argument is the constant coefficient of the
     * polynomial.
     *
     * @param root0 If the polynomial has real roots, this reference
     * argument is used to return the first root.
     *
     * @param root1 If the polynomial has real roots, this reference
     * argument is used to return the second root.
     *
     * @param valid If the polynomial has real roots, this reference
     * argument is set to true.  If the polynomial does not have real
     * roots, this polynomial is set to false.
     *
     * @param checkValidity The first thing solveQuadratic does is check
     * to see if the polynomial has real roots.  By setting
     * checkValidity to false, the user can disable this check.  It
     * makes sense to do this only if the polynomial is already known to
     * have real roots, and execution time is at a premium.
     */
    template <class Type>
    void
    solveQuadratic(Type c0, Type c1, Type c2,
                   Type& root0, Type& root1, bool& valid,
                   bool checkValidity=true);
#endif


    /**
     * This function computes the standard deviation of the elements
     * of its argument, and allows the user to specify the precision
     * with which the computation is carried out.  The computation is
     * performed using the type specified by the second template
     * argument.  This routine eventually calls the function
     * mean(Array1D const&, ...).  This is a little inefficient, since
     * you are probably already calling mean(...).  If this redundant
     * computation is not OK, then we should add a class called
     * "SignalMoments" (or something similar) which will compute both
     * mean and standard deviation without any wasted cycles.  Perhaps
     * by the time you read this, such a class will have already been
     * implemented.
     *
     * @param array0 This argument is the input array, from the elements
     * of which the standard deviation will be computed.
     *
     * @return The return value is the mean of the elements of array0.
     */
    template <class Type0, class Type1>
    inline Type1
    standardDeviation(Array1D<Type0> const& array0);


    /**
     * This function computes the sum of the elements of its argument.
     * The summation is accumulated into a variable of type Type2,
     * allowing the user to control the precision of the internal
     * summation.
     *
     * @param array0 This argument is the array to be summed.
     *
     * @return The summation of all the elements of array0.
     */
    template <class Type2, class Type>
    Type2
    sum(Array1D<Type> const& array0);


    /**
     * This function computes the sum of those elements of its
     * argument which lie within a rectangular region of interest.
     * The summation is accumulated into a variable of type Type2,
     * allowing the user to control the precision of the internal
     * summation.
     *
     * @param array0 This argument is the array to be summed.
     *
     * @param upperLeftCorner This argument specifies the upper left
     * corner of the rectangular region to be summed.  The summed
     * region will include the array element corresponding to
     * upperLeftCorner.
     *
     * @param lowerRightCorner This argument specifies the lower right
     * corner of the rectangular region to be summed.  The summed
     * region will stop one row/column short of lowerRightCorner.
     *
     * @return The summation of all the elements in the region of
     * interest.
     */
    template <class Type2, class Type>
    Type2
    sum(Array2D<Type> const& array0,
        Index2D const& upperLeftCorner,
        Index2D const& lowerRightCorner);


    /**
     * This function returns an array made up of only those elements
     * of dataArray that correspond to indices in indexArray.  For
     * example, in the code below, resultArray should end up with
     * four elements ([2.0, 1.0, 1.0, 4.0]):
     *
     * @code
     *   Array1D<double> dataArray("[0.0, 1.0, 2.0, 3.0, 4.0, 5.0]");
     *   Array1D<unsigned int> indexArray("[2, 1, 1, 4]");
     *   Array1D<double> resultArray = take(dataArray, indexArray);
     * @endCode
     *
     * It is an error if indexArray contains indices that are not
     * valid for dataArray.
     *
     * @param dataArray This argument is the array from which to draw
     * elements.
     *
     * @param indexArray This argument contains, in order, the indices
     * of the elements that should be included in the output.
     *
     * @return The return value is an array containing those elements
     * selected by indexArray.
     */
    template <class Type, class IntegralType>
    Array1D<Type>
    take(Array1D<Type> const& dataArray,
         Array1D<IntegralType> const& indexArray);


    /**
     * This function returns an array made up of only those elements
     * of dataArray that correspond to indices in indexArray.  Its
     * operation is just like that of take(Array1D const&, Array1D
     * const&), with the additional detail that the input array is
     * flattened before the take operation.  For example, in the code
     * below, resultArray should end up with four elements ([2.0, 1.0,
     * 1.0, 4.0]):
     *
     * @code
     *   Array2D<double> dataArray("[[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]]");
     *   Array1D<unsigned int> indexArray("[2, 1, 1, 4]");
     *   Array1D<double> resultArray = take(dataArray, indexArray);
     * @endCode
     *
     * It is an error if indexArray contains indices that are less
     * than zero or greater than (dataArray.size() - 1).
     *
     * @param dataArray This argument is the array from which to draw
     * elements.
     *
     * @param indexArray This argument contains, in order, the indices
     * of the elements that should be included in the output.
     *
     * @return The return value is an array containing those elements
     * selected by indexArray.
     */
    template <class Type, class IntegralType>
    inline Array1D<Type>
    take(Array2D<Type> const& dataArray,
         Array1D<IntegralType> const& indexArray);


    /**
     * This function works just like take(Array2D const&, Array1D
     * const&), with the exception that the input array is not
     * flattened, and entire rows (or columns) are selected.  For
     * example, in the code below, resultArray0 should end up with two
     * rows and three columns ([[0.0, 1.0, 2.0], [6.0, 7.0, 8.0]]),
     * while resultArray1 should end up with four rows and two columns
     * ([[0.0, 2.0], [3.0, 5.0], [6.0, 8.0], [9.0, 11.0]]):
     *
     * @code
     *   Array2D<double> dataArray("[[0.0, 1.0, 2.0], "
     *                             " [3.0, 4.0, 5.0], "
     *                             " [6.0, 7.0, 8.0], "
     *                             " [9.0, 10.0, 11.0]]");
     *   Array1D<unsigned int> indexArray("[0, 2]");
     *   Array2D<double> resultArray0 = take(dataArray, indexArray, 0);
     *   Array2D<double> resultArray1 = take(dataArray, indexArray, 1);
     * @endCode
     *
     * It is an error if axis == 0 and indexArray contains indices
     * that are less than zero or greater than (dataArray.rows() - 1),
     * or if axis == 1 and indexArray contains indices that are less
     * than zero or greater than (dataArray.columns() - 1).
     *
     * @param dataArray This argument is the array from which to draw
     * elements.
     *
     * @param indexArray This argument contains, in order, the indices
     * of the elements that should be included in the output.
     *
     * @param axis This argument specifies the axis over which to
     * select values.  If axis == 0, then individual rows will be
     * taken.  If axis == 1, then individual columns will be taken.
     * Behavior is undefined if axis is not equal to either 0 or 1.
     *
     * @return The return value is an array containing those elements
     * selected by indexArray.
     */
    template <class Type, class IntegralType>
    Array2D<Type>
    take(Array2D<Type> const& dataArray,
         Array1D<IntegralType> const& indexArray,
         unsigned int axis);


    /**
     * This function computes the variance, of the elements of its
     * argument, and allows the user to specify the precision with which
     * the computation is carried out.  The computation is performed
     * using the type specified by the second argument to the variance()
     * function call.  This routine eventually calls the function
     * mean(Array1D const&, ...).  This is a little inefficient, since
     * you are probably already calling mean(...).  If this redundant
     * computation is not OK, then we should add a class called
     * "SignalMoments" (or something similar) which will compute both
     * mean and variance without any wasted cycles.  Perhaps by the time
     * you read this, such a class will have already been implemented.
     *
     * @param array0 This argument is the input array, from the elements
     * of which the variance will be computed.
     *
     * @param typeTag This argument specifies the type which
     * will be used to perform the computation, and also specifies the
     * return type of the function.
     *
     * @return The return value is the mean of the elements of array0.
     */
    template <class Type0, class Type1>
    inline Type1
    variance(Array1D<Type0> const& array0);


    /**
     * This function returns an Array1D of the specified size and type
     * in which the value of every element is zero.
     *
     * @param size This argument specifies the number of elements in the
     * returned array.
     *
     * @return An Array1D instance in which each element is set to zero.
     */
    template <class Type>
    inline Array1D<Type>
    zeros(size_t size);


    /**
     * This function returns an Array2D of the specified size and type
     * in which the value of every element is zero.
     *
     * @param rows This argument specifies the number of rows in the
     * returned array.
     *
     * @param columns This argument specifies the number of columns in the
     * returned array.
     *
     * @return An Array2D instance in which each element is set to zero.
     */
    template <class Type>
    inline Array2D<Type>
    zeros(size_t rows, size_t columns);


    /**
     * This function returns an Array3D of the specified size and type
     * in which the value of every element is zero.
     *
     * @param shape0 This argument specifies the number of slices in the
     * returned array.
     *
     * @param shape1 This argument specifies the number of rows in the
     * returned array.
     *
     * @param shape2 This argument specifies the number of columns in the
     * returned array.
     *
     * @return An Array3D instance in which each element is set to zero.
     */
    template <class Type>
    inline Array3D<Type>
    zeros(size_t shape0, size_t shape1, size_t shape2);

  } // namespace numeric

} // namespace brick

#include <brick/numeric/utilities_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_UTILITIES_HH */
