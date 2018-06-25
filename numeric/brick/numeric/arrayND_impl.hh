/**
***************************************************************************
* @file brick/numeric/arrayND_impl.hh
*
* Header file defining ArrayND class template.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#ifndef BRICK_NUMERIC_ARRAYND_IMPL_HH
#define BRICK_NUMERIC_ARRAYND_IMPL_HH

// This file is included by arrayND.hh, and should not be directly included
// by user code, so no need to include arrayND.hh here.
//
// #include <brick/numeric/arrayND.hh>

#include <algorithm>
#include <numeric>
#include <sstream>
#include <vector>
#include <brick/numeric/numericTraits.hh>

namespace brick {

  namespace numeric {

    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>::
    ArrayND()
      : m_shape(),
        m_storage(),
        m_strideArray()
    {
      // Empty.
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>::
    ArrayND(size_t dimension, size_t const shape[])
      : m_shape(dimension),
        m_storage(),
        m_strideArray()
    {
      m_shape.copy(shape);
      m_storage.reinit(this->computeSize(m_shape));
    }


    // Construct from an initialization string.
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>::
    ArrayND(const Array1D<size_t>& shape)
      : m_shape(shape.copy()),
        m_storage(this->computeSize(m_shape)),
        m_strideArray(this->computeStride(m_shape))
    {
      // Empty.
    }


    /* When copying from a ArrayND do a shallow copy */
    /* Update reference count if the array we're copying has */
    /* valid data. */
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>::
    ArrayND(const ArrayND<Dimension, Type>& source)
      : m_shape(source.m_shape.copy()),            // Deep copy.
        m_storage(source.m_storage),               // Shallow copy.
        m_strideArray(source.m_strideArray.copy()) // Deep copy.
    {
      // Empty.
    }


    // Construct from an initialization string.
    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>::
    ArrayND(const std::string& inputString)
    {
      std::istringstream stream(inputString);
      stream >> *this;
      if(!stream) {
        BRICK_THROW(brick::common::ValueException,
                    "ArrayND::ArrayND(std::string const&)",
                    "Malformed input string.");
      }
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>::
    ~ArrayND()
    {
      // Empty.
    }


    template <size_t Dimension, class Type>
    inline void ArrayND<Dimension, Type>::
    checkDimension(Array1D<size_t> const& shape) const
    {
      // Hmm.  Seems like this should happen even if we're optimizing.
// #ifdef BRICK_NUMERIC_CHECKBOUNDS
      if(shape.size() != m_shape.size()) {
        std::ostringstream message;
        message << "Dimensionality mismatch: required dimensionality is "
                << shape.size()
                << " while *this has dimension " << m_shape.size() << ".";
        BRICK_THROW(brick::common::IndexException, "ArrayND::checkDimension()",
                    message.str().c_str());
      }
      for(size_t dimension = 0; dimension < m_shape.size(); ++dimension) {
        if(m_shape[dimension] != shape[dimension]) {
          std::ostringstream message;
          message << "Shape mismatch: required shape is " << shape
                << " while *this has shape " << m_shape << ".";
          BRICK_THROW(brick::common::IndexException,
                      "ArrayND::checkDimension()",
                      message.str().c_str());
        }
      }
// #endif
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> ArrayND<Dimension, Type>::
    copy() const
    {
      ArrayND<Dimension, Type> newArray(m_shape);
      newArray.copy(*this);
      return newArray;
    }


    template <size_t Dimension, class Type> template <class Type2>
    void ArrayND<Dimension, Type>::
    copy(const ArrayND<Dimension, Type2>& source)
    {
      this->checkDimension(source.getShape());
      if(m_storage.size() != 0) {
        this->copy(source.data());
      }
    }


    template <size_t Dimension, class Type> template <class Type2>
    void
    ArrayND<Dimension, Type>::
    copy(const Type2* dataPtr)
    {
      m_storage.copy(dataPtr);
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>&
    ArrayND<Dimension, Type>::
    operator=(Type val)
    {
      m_storage = val;
      return *this;
    }


    template <size_t Dimension, class Type>
    size_t
    ArrayND<Dimension, Type>::
    flattenIndex(const Array1D<size_t>& indexArray) const
    {
      brick::numeric::Array1D<size_t>::const_iterator endMinus = indexArray.end();
      --endMinus;
      return std::inner_product(indexArray.begin(), endMinus,
                                m_strideArray.begin(), *endMinus);
    }


    template <size_t Dimension, class Type>
    void ArrayND<Dimension, Type>::
    reinit(Array1D<size_t> const& shape)
    {
      m_shape = shape.copy();
      m_storage.reinit(this->computeSize(m_shape));
      m_strideArray = this->computeStride(m_shape);
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>& ArrayND<Dimension, Type>::
    operator=(const ArrayND<Dimension, Type>& source)
    {
      // Check for self-assignment
      if(&source != this) {
        m_shape = source.m_shape.copy();             // Deep copy.
        m_storage = source.m_storage;                // Shallow copy.
        m_strideArray = source.m_strideArray.copy(); // Deep copy.
      }
      return *this;
    }


    template <size_t Dimension, class Type> template <class Type2>
    ArrayND<Dimension, Type>&
    ArrayND<Dimension, Type>::
    operator+=(const ArrayND<Dimension, Type2>& arg)
    {
      // Let m_storage worry about array size errors.
      m_storage += arg.m_storage;
      return *this;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>&
    ArrayND<Dimension, Type>::
    operator+=(const Type arg)
    {
      m_storage += arg;
      return *this;
    }


    template <size_t Dimension, class Type> template <class Type2>
    ArrayND<Dimension, Type>&
    ArrayND<Dimension, Type>::
    operator-=(const ArrayND<Dimension, Type2>& arg)
    {
      // Let m_storage worry about array size errors.
      m_storage -= arg.m_storage;
      return *this;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>&
    ArrayND<Dimension, Type>::
    operator-=(const Type arg)
    {
      m_storage -= arg;
      return *this;
    }


    template <size_t Dimension, class Type> template <class Type2>
    ArrayND<Dimension, Type>&
    ArrayND<Dimension, Type>::
    operator*=(const ArrayND<Dimension, Type2>& arg)
    {
      // Let m_storage worry about array size errors.
      m_storage *= arg.m_storage;
      return *this;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>&
    ArrayND<Dimension, Type>::
    operator*=(const Type arg)
    {
      m_storage *= arg;
      return *this;
    }


    template <size_t Dimension, class Type> template <class Type2>
    ArrayND<Dimension, Type>&
    ArrayND<Dimension, Type>::
    operator/=(const ArrayND<Dimension, Type2>& arg)
    {
      // Let m_storage worry about array size errors.
      m_storage /= arg.m_storage;
      return *this;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type>&
    ArrayND<Dimension, Type>::
    operator/=(const Type arg)
    {
      m_storage /= arg;
      return *this;
    }


    template <size_t Dimension, class Type>
    size_t
    ArrayND<Dimension, Type>::
    computeSize(Array1D<size_t> const& shape) const
    {
      return std::accumulate(shape.begin(), shape.end(), size_t(1),
                             std::multiplies<size_t>());
    }


    template <size_t Dimension, class Type>
    Array1D<size_t>
    ArrayND<Dimension, Type>::
    computeStride(Array1D<size_t> const& shape) const
    {
      Array1D<size_t> strideArray(shape.size());
      strideArray[shape.size() - 1] = 1;
      for(size_t ii = shape.size() - 1; ii != 0; --ii) {
        strideArray[ii - 1] = shape[ii] * strideArray[ii];
      }
      return strideArray;
    }


    /* ========== Non-member functions ========== */

    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator+(const ArrayND<Dimension, Type>& array0,
                                       const ArrayND<Dimension, Type>& arrayN)
    {
      if(array0.size() != arrayN.size()) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 has " << array0.size()
                << " elements, while arrayN has " << arrayN.size()
                << " elements.";
        BRICK_THROW(brick::common::ValueException, "ArrayND::operator+()",
                    message.str().c_str());
      }
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), arrayN.begin(),
                     result.begin(), std::plus<Type>());
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator-(const ArrayND<Dimension, Type>& array0,
                            const ArrayND<Dimension, Type>& arrayN)
    {
      if(array0.size() != arrayN.size()) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 has " << array0.size()
                << " elements, while arrayN has " << arrayN.size()
                << " elements.";
        BRICK_THROW(brick::common::ValueException, "ArrayND::operator-()",
                    message.str().c_str());
      }
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), arrayN.begin(),
                     result.begin(), std::minus<Type>());
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator*(const ArrayND<Dimension, Type>& array0,
                            const ArrayND<Dimension, Type>& arrayN)
    {
      if(array0.size() != arrayN.size()) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 has " << array0.size()
                << " elements, while arrayN has " << arrayN.size()
                << " elements.";
        BRICK_THROW(brick::common::ValueException,
                    "ArrayND::operator*()", message.str().c_str());
      }
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), arrayN.begin(),
                     result.begin(), std::multiplies<Type>());
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator/(const ArrayND<Dimension, Type>& array0,
                            const ArrayND<Dimension, Type>& arrayN)
    {
      if(array0.size() != arrayN.size()) {
        std::ostringstream message;
        message << "Array sizes do not match.  Array0 has " << array0.size()
                << " elements, while arrayN has " << arrayN.size()
                << " elements.";
        BRICK_THROW(brick::common::ValueException,
                    "ArrayND::operator/()", message.str().c_str());
      }
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), arrayN.begin(),
                     result.begin(), std::divides<Type>());
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator+(const ArrayND<Dimension, Type>& array0, Type scalar)
    {
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::plus<Type>(), scalar));
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator-(const ArrayND<Dimension, Type>& array0, Type scalar)
    {
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::minus<Type>(), scalar));
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator*(const ArrayND<Dimension, Type>& array0, Type scalar)
    {
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::multiplies<Type>(), scalar));
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator/(const ArrayND<Dimension, Type>& array0, Type scalar)
    {
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::divides<Type>(), scalar));
      return result;
    }


    template <size_t Dimension, class Type>
    inline ArrayND<Dimension, Type> operator+(Type scalar, const ArrayND<Dimension, Type>& array0)
    {
      return array0 + scalar;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator-(Type scalar, const ArrayND<Dimension, Type>& array0)
    {
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind1st(std::minus<Type>(), scalar));
      return result;
    }


    template <size_t Dimension, class Type>
    inline ArrayND<Dimension, Type> operator*(Type scalar, const ArrayND<Dimension, Type>& array0)
    {
      return array0 * scalar;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, Type> operator/(Type scalar, const ArrayND<Dimension, Type>& array0)
    {
      ArrayND<Dimension, Type> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind1st(std::divides<Type>(), scalar));
      return result;
    }


    // Elementwise comparison of an ArrayND with a constant.
    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator==(const ArrayND<Dimension, Type>& array0, const Type arg)
    {
      ArrayND<Dimension, bool> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.data(),
                     std::bind2nd(std::equal_to<Type>(), arg));
      return result;
    }


    // Elementwise comparison of an ArrayND with another array.
    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator==(const ArrayND<Dimension, Type>& array0,
               const ArrayND<Dimension, Type>& arrayN)
    {
      array0.checkDimension(arrayN.getShape());
      ArrayND<Dimension, bool> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), arrayN.begin(),
                     result.begin(), std::equal_to<Type>());
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator>(const ArrayND<Dimension, Type>& array0, const Type arg)
    {
      ArrayND<Dimension, bool> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::greater<Type>(), arg));
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator>=(const ArrayND<Dimension, Type>& array0, const Type arg)
    {
      ArrayND<Dimension, bool> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::greater_equal<Type>(), arg));
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator<(const ArrayND<Dimension, Type>& array0, const Type arg)
    {
      ArrayND<Dimension, bool> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::less<Type>(), arg));
      return result;
    }


    template <size_t Dimension, class Type>
    ArrayND<Dimension, bool>
    operator<=(const ArrayND<Dimension, Type>& array0, const Type arg)
    {
      ArrayND<Dimension, bool> result(array0.getShape());
      std::transform(array0.begin(), array0.end(), result.begin(),
                     std::bind2nd(std::less_equal<Type>(), arg));
      return result;
    }


    template <size_t Dimension, class Type>
    std::ostream& operator<<(std::ostream& stream,
                             const ArrayND<Dimension, Type>& array0)
    {
      if (!stream){
        return stream;
      }
      stream << "ArrayND {\n"
             << "  shape: " << array0.getShape() << "\n"
             << "  data: " << array0.ravel() << "\n"
             << "}";
      return stream;
    }


    // Sets the value of an ArrayND instance from a std::istream.
    template <size_t Dimension, class Type>
    std::istream& operator>>(std::istream& stream,
                             ArrayND<Dimension, Type> & array0)
    {
      if (!stream){
        return stream;
      }

      common::Expect::FormatFlag flags = common::Expect::SkipWhitespace();
      stream >> common::Expect("ArrayND", flags)
             >> common::Expect("{", flags);

      Array1D<size_t> shape;
      stream >> common::Expect("shape:", flags);
      stream >> shape;

      Array1D<Type> data;
      stream >> common::Expect("data:", flags);
      stream >> data;

      if(stream) {
        array0.reinit(shape);
        array0.copy(data.data());
      }
      return stream;
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifdef BRICK_NUMERIC_ARRAYND_IMPL_HH */
