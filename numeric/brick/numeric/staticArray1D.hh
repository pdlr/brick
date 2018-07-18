/**
***************************************************************************
* @file brick/numeric/staticArray1D.hh
*
* Header file declaring StaticArray1D class.
*
* Copyright (C) 2001-2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_STATICARRAY1D_HH
#define BRICK_NUMERIC_STATICARRAY1D_HH

#include <iostream>
#include <string>
#include <brick/common/exception.hh>

namespace brick {

  namespace numeric {

    /**
     ** The StaticArray1D class template represents a 1D array of
     ** arbitrary type.  All memory used by this class is obtained from
     ** the stack.  Unlike Array1D, StaticArray1D has fixed size
     ** (determined by a template parameter at compile time) and has
     ** deep copy semantics, so that if you write:
     **
     **       staticArray1 = staticArray2;
     **
     ** then the data from staticArray2 will be deep copied into staticArray1.
     **/
    template <class Type, size_t Size>
    class StaticArray1D {
    public:
      /* ======== Public typedefs ======== */

      /**
       ** Typedef for value_type describes the contents of the staticArray.
       **/
      typedef Type value_type;


      /**
       ** Typedef for iterator type helps with standard library interface.
       **/
      typedef Type* iterator;


      /**
       ** Typedef for const_iterator type helps with standard library
       ** interface.
       **/
      typedef const Type* const_iterator;


      /* ======== Public member functions ======== */

      /**
       * Default constructor.
       */
      StaticArray1D();


      /**
       * Constructs an StaticArray1D by parsing an input string.  For example,
       * the code line
       *
       *   StaticArray1D<double, 3> testStaticArray("[1.0, 2.0, 3.0]");
       *
       * would construct a 3 element staticArray.
       *
       * @param inputString A formatted string which could be reasonably
       * expected to parse into the desired staticArray values.
       */
      explicit
      StaticArray1D(const std::string& inputString);


      /**
       * The copy constructor does a deep copy.
       *
       * @param source The StaticArray1D instance to be copied.
       */
      StaticArray1D(const StaticArray1D<Type, Size> &source);


      /**
       * Destroys the StaticArray1D instance.
       */
      ~StaticArray1D();


      /**
       * Return begin() iterator for Standard Library algorithms.
       *
       * @return Iterator pointing to the first element of the staticArray.
       */
      iterator
      begin() {return &(m_dataArray[0]);}


      /**
       * Return begin() const_iterator for Standard Library algorithms.
       *
       * @return Const iterator pointing to the first element of the
       * staticArray.
       */
      const_iterator
      begin() const {return &(m_dataArray[0]);}


      /**
       * Return end() iterator for Standard Library algorithms.
       *
       * @return Iterator pointing just past the last element of
       * the staticArray.
       */
      iterator
      end() {this->begin() + Size;}


      /**
       * Return end() const_iterator for Standard Library algorithms.
       *
       * @return Const iterator pointing just past the last element of
       * the staticArray.
       */
      const_iterator
      end() const {this->begin() + Size;}


      /**
       * Returns the number of elements in the staticArray.  This is a synonym
       * for size().
       *
       * @return The number of elements in the staticArray.
       */
      size_t
      length() const {return this->size();}


      /**
       * Returns the number of elements in the staticArray.  This is a synonym
       * for length().
       *
       * @return The number of elements in the staticArray.
       */
      size_t
      size() const {return Size;}


      /**
       * Assignment operator deep copies the contents of source.
       *
       * @param source The StaticArray1D instance to be copied.
       *
       * @return Reference to *this.
       */
      StaticArray1D<Type, Size>&
      operator=(const StaticArray1D<Type, Size>& source);


      /**
       * Assign value to every element in the staticArray.
       *
       * @param value The value to be copied.
       *
       * @return Reference to *this.
       */
      StaticArray1D<Type, Size>&
      operator=(Type value);


      /**
       * Returns the (index)th element of the staticArray by reference.
       *
       * @param index Indicates the selected element.
       *
       * @return Reference to the (index)th element of the staticArray.
       */
      inline Type&
      operator()(size_t index) {
        this->checkBounds(index);
        return m_dataArray[index];
      }


      /**
       * Returns the (index)th element of the staticArray by value.
       *
       * @param index Indicates the selected element.
       *
       * @return Value of the (index)th element of the staticArray.
       */
      Type operator()(size_t index) const {
        this->checkBounds(index);
        return m_dataArray[index];
      }


      /**
       * Returns the (index)th element of the staticArray by reference.
       * Synonymous with operator()().
       *
       * @param index Indicates the selected element.
       *
       * @return Reference to the (index)th element of the staticArray.
       */
      Type& operator[](size_t index) {return this->operator()(index);}


      /**
       * Returns the (index)th element of the staticArray by value.
       * Synonymous with operator()() const.
       *
       * @param index Indicates the selected element.
       *
       * @return Value of the (index)th element of the staticArray.
       */
      Type operator[](size_t index) const {return this->operator()(index);}


      /**
       * Increments each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg StaticArray1D of values to be added to the elements of
       * *this.
       *
       * @exception ValueException on incompatible staticArray sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      StaticArray1D<Type, Size>&
      operator+=(const StaticArray1D<Type2, Size>& arg);


      /**
       * Increments each element of *this by a constant.
       *
       * @param arg Value by which staticArray elements will be incremented.
       *
       * @return Reference to *this.
       */
      StaticArray1D<Type, Size>&
      operator+=(const Type arg);


      /**
       * Decrements each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg StaticArray1D of values to be subtracted from the elements
       * of *this.
       *
       * @exception ValueException on incompatible staticArray sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      StaticArray1D<Type, Size>&
      operator-=(const StaticArray1D<Type2, Size>& arg);


      /**
       * Decrements each element of *this by a constant.
       *
       * @param arg Value by which staticArray elements will be decremented.
       *
       * @return Reference to *this.
       */
      StaticArray1D<Type, Size>&
      operator-=(const Type arg);


      /**
       * Multiplies each element of *this by the value of the
       * corresponding element of arg.
       *
       * @param arg StaticArray1D of values by which the elements of *this
       * should be multiplied.
       *
       * @exception ValueException on incompatible staticArray sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      StaticArray1D<Type, Size>&
      operator*=(const StaticArray1D<Type2, Size>& arg);


      /**
       * Multiplies each element of *this by a constant.
       *
       * @param arg Value by which staticArray elements will be multiplied.
       *
       * @return Reference to *this.
       */
      StaticArray1D<Type, Size>&
      operator*=(const Type arg);


      /**
       * Divides each element of *this by the value of the corresponding
       * element of arg.
       *
       * @param arg StaticArray1D of values by which the elements of *this
       * should be divided.
       *
       * @exception ValueException on incompatible staticArray sizes
       *
       * @return Reference to *this.
       */
      template <class Type2>
      StaticArray1D<Type, Size>&
      operator/=(const StaticArray1D<Type2, Size>& arg);


      /**
       * Divides each element of *this by a constant.
       *
       * @param arg Value by which staticArray elements will be divided.
       *
       * @return Reference to *this.
       */
      StaticArray1D<Type, Size>&
      operator/=(const Type arg);


      /**
       * Elementwise comparison with a constant.
       *
       * @param arg Value to which staticArray elements should be compared.
       *
       * @return StaticArray1D<bool, Size> in which each element has value "true" if
       * the corresponding element of *this is greater than arg.
       */
      StaticArray1D<bool, Size>
      operator>(const Type arg) const;


      /**
       * Elementwise comparison with a constant.
       *
       * @param arg Value to which staticArray elements should be compared.
       *
       * @return StaticArray1D<bool, Size> in which each element has value "true" if
       * the corresponding element of *this is greater than or equal to arg.
       */
      StaticArray1D<bool, Size>
      operator>=(const Type arg) const;


      /**
       * Elementwise comparison with a constant.
       *
       * @param arg Value to which staticArray elements should be compared.
       *
       * @return StaticArray1D<bool, Size> in which each element has value "true" if
       * the corresponding element of *this is less than arg.
       */
      StaticArray1D<bool, Size>
      operator<(const Type arg) const;


      /**
       * Elementwise comparison with a constant.
       *
       * @param arg Value to which staticArray elements should be compared.
       *
       * @return StaticArray1D<bool, Size> in which each element has value "true" if
       * the corresponding element of *this is less than or equal to arg.
       */
      StaticArray1D<bool, Size>
      operator<=(const Type arg) const;

    private:
      /* ======== Private member functions ======== */

      /**
       * Optionally throw an exception if index is beyond the range of
       * this staticArray.
       *
       * @param index The index to check.
       */
      inline void
      checkBounds(size_t index) const;


      // Constants to help with formatting.  We use the initialization
      // on first use paradigm for the string constants to avoid
      // headaches.

      /**
       * Static constant describing how the string representation of an
       * StaticArray1D should start.
       */
      static const std::string& ioIntro();


      /**
       * Static constant describing how the string representation of an
       * StaticArray1D should end.
       */
      static const std::string& ioOutro();


      /**
       * Static constant describing how the the data portion of the
       * string representation of an StaticArray1D should start.
       */
      static const char ioOpening = '[';


      /**
       * Static constant describing how the the data portion of the
       * string representation of an StaticArray1D should end.
       */
      static const char ioClosing = ']';


      /**
       * Static constant describing how individual elements should be
       * separated in the string representation of StaticArray1D.
       */
      static const char ioSeparator = ',';


      /* ======== Private data members ======== */
      Type m_dataArray[Size];

    };


    /* ================= Non-member functions ================= */

    /**
     * Elementwise addition of StaticArray1D instances.
     *
     * @param staticArray0 First argument for addition.
     *
     * @param staticArray1 Second argument for addition.
     *
     * @exception ValueException on incompatible staticArray sizes
     *
     * @return StaticArray1D instance in which the value of each element is
     * the sum of the values of the corresponding elements of the two
     * StaticArray1D arguments.
     */
    template <class Type>
    StaticArray1D<Type, Size>
    operator+(const StaticArray1D<Type, Size>& staticArray0,
              const StaticArray1D<Type, Size>& staticArray1);


    /**
     * Elementwise subtraction of StaticArray1D instances.
     *
     * @param staticArray0 First argument for subtraction.
     *
     * @param staticArray1 Second argument for subtraction.
     *
     * @exception ValueException on incompatible staticArray sizes
     *
     * @return StaticArray1D instance in which the value of each element is
     * the difference of the values of the corresponding elements of the two
     * StaticArray1D arguments.
     */
    template <class Type>
    StaticArray1D<Type, Size>
    operator-(const StaticArray1D<Type, Size>& staticArray0,
              const StaticArray1D<Type, Size>& staticArray1);


    /**
     * Elementwise multiplication of StaticArray1D instances.
     *
     * @param staticArray0 First argument for multiplication.
     *
     * @param staticArray1 Second argument for multiplication.
     *
     * @exception ValueException on incompatible staticArray sizes
     *
     * @return StaticArray1D instance in which the value of each element is
     * the product of the values of the corresponding elements of the two
     * StaticArray1D arguments.
     */
    template <class Type>
    StaticArray1D<Type, Size>
    operator*(const StaticArray1D<Type, Size>& staticArray0,
              const StaticArray1D<Type, Size>& staticArray1);


    /**
     * Elementwise division of StaticArray1D instances.
     *
     * @param staticArray0 First argument for division.
     *
     * @param staticArray1 Second argument for division.
     *
     * @exception ValueException on incompatible staticArray sizes
     *
     * @return StaticArray1D instance in which the value of each element is
     * the dividend of the values of the corresponding elements of the two
     * StaticArray1D arguments.
     */
    template <class Type>
    StaticArray1D<Type, Size>
    operator/(const StaticArray1D<Type, Size>& staticArray0,
              const StaticArray1D<Type, Size>& staticArray1);


    /**
     * Addition of StaticArray1D and scalar.
     *
     * @param staticArray StaticArray1D argument of the addition.
     *
     * @param scalar Scalar argument of the addition.
     *
     * @return StaticArray1D instance in which the value of each element is
     * the sum of the corresponding element of the StaticArray1D argument and
     * the scalar argument.
     */
    template <class Type>
    StaticArray1D<Type, Size>
    operator+(const StaticArray1D<Type, Size>& staticArray, Type scalar);


    /**
     * Subtraction of StaticArray1D and scalar.
     *
     * @param staticArray0 StaticArray1D argument of the subtraction.
     *
     * @param scalar Scalar argument of the subtraction.
     *
     * @return StaticArray1D instance in which the value of each element is
     * the difference of the corresponding element of the StaticArray1D
     * argument and the scalar argument.
     */
    template <class Type>
    StaticArray1D<Type, Size>
    operator-(const StaticArray1D<Type, Size>& staticArray0, Type scalar);


    /**
     * Multiplication of StaticArray1D and scalar.
     *
     * @param staticArray0 StaticArray1D argument of the multiplication.
     *
     * @param scalar Scalar argument of the multiplication.
     *
     * @return StaticArray1D instance in which the value of each element is
     * the product of the corresponding element of the StaticArray1D argument
     * and the scalar argument.
     */
    template <class Type>
    StaticArray1D<Type, Size>
    operator*(const StaticArray1D<Type, Size>& staticArray0, Type scalar);


    /**
     * Division of StaticArray1D and scalar.
     *
     * @param staticArray0 StaticArray1D argument of the division.
     *
     * @param scalar Scalar argument of the division.
     *
     * @return StaticArray1D instance in which the value of each element is
     * the difference of the corresponding element of the StaticArray1D
     * argument and the scalar argument.
     */
    template <class Type>
    StaticArray1D<Type, Size>
    operator/(const StaticArray1D<Type, Size>& staticArray0, Type scalar);


    /**
     * Addition of scalar and StaticArray1D.
     *
     * @param scalar Scalar argument of the addition.
     *
     * @param staticArray0 StaticArray1D argument of the addition.
     *
     * @return StaticArray1D instance in which the value of each element is
     * the sum of the scalar argument and the corresponding element of
     * the StaticArray1D argument.
     */
    template <class Type>
    inline StaticArray1D<Type, Size>
    operator+(Type scalar, const StaticArray1D<Type, Size>& staticArray0);


    /**
     * Subtraction of scalar and StaticArray1D.
     *
     * @param scalar Scalar argument of the subtraction.
     *
     * @param staticArray0 StaticArray1D argument of the subtraction.
     *
     * @return StaticArray1D instance in which the value of each element is
     * the difference of the scalar argument and the corresponding
     * element of the StaticArray1D argument.
     */
    template <class Type>
    inline StaticArray1D<Type, Size>
    operator-(Type scalar, const StaticArray1D<Type, Size>& staticArray0);


    /**
     * Multiplication of scalar and StaticArray1D.
     *
     * @param scalar Scalar argument of the multiplication.
     *
     * @param staticArray0 StaticArray1D argument of the multiplication.
     *
     * @return StaticArray1D instance in which the value of each element is
     * the product of the scalar argument and the corresponding element
     * of the StaticArray1D argument.
     */
    template <class Type>
    inline StaticArray1D<Type, Size>
    operator*(Type scalar, const StaticArray1D<Type, Size>& staticArray0);


    /**
     * Division of scalar and StaticArray1D.
     *
     * @param scalar Scalar argument of the division.
     *
     * @param staticArray0 StaticArray1D argument of the division.
     *
     * @return StaticArray1D instance in which the value of each element is
     * the dividend of the scalar argument and the corresponding element
     * of the StaticArray1D argument.
     */
    template <class Type>
    inline StaticArray1D<Type, Size>
    operator/(Type scalar, const StaticArray1D<Type, Size>& staticArray0);


    /**
     * Outputs a text representation of an StaticArray1D instance to a
     * std::ostream.  The output format looks like this:
     *
     * StaticArray1D([1, 2, 4, 8, 16])
     *
     * Where the staticArray elements are output using
     * operator<<(std::ostream&, const Type&)
     * and each element is separated from its neighbors by a comma and
     * whitespace.
     *
     * @param stream Reference to the the output stream.
     *
     * @param staticArray0 const Reference to the StaticArray1D to be output.
     *
     * @exception IOException on invalid stream.
     *
     * @return Reference to output stream.
     */
    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream,
               const StaticArray1D<Type, Size>& staticArray0);


    /**
     * Sets the value of an StaticArray1D instance from a std::istream.
     * The input format is as described for
     * operator<<(std::ostream&, const StaticArray1D<Type, Size>&) above.
     *
     * @param stream Reference to the the input stream.
     *
     * @param staticArray0 Reference to the StaticArray1D which will
     * take the input.
     *
     * @return Reference to input stream.
     */
    template <class Type>
    std::istream&
    operator>>(std::istream& stream,
               StaticArray1D<Type, Size>& staticArray0);

  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/staticArray1D_impl.hh>

#endif /* #ifdef BRICK_NUMERIC_STATICARRAY1D_HH */
