/**
***************************************************************************
* @file brick/numeric/vector3D_impl.hh
*
* Header file defining inline and template functions declared in
* vector3D.hh.
*
* Copyright (C) 2000-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_VECTOR3D_IMPL_HH
#define BRICK_NUMERIC_VECTOR3D_IMPL_HH

// This file is included by vector3D.hh, and should not be directly included
// by user code, so no need to include vector3D.hh here.
//
// #include <brick/numeric/vector3D.hh>

namespace brick {

  namespace numeric {

    //  Default constructor initializes to (0, 0, 0).
    template <class Type>
    Vector3D<Type>::
    Vector3D()
      : m_x(0), m_y(0), m_z(0)
    {
      // Empty.
    }


    // Explicitly sets 3D coordinates.
    template <class Type>
    Vector3D<Type>::
    Vector3D(Type const& xCoord, Type const& yCoord, Type const& zCoord)
      : m_x(xCoord), m_y(yCoord), m_z(zCoord) {}


    // Explicitly sets 3D homogeneous coordinates.
    template <class Type>
    Vector3D<Type>::
    Vector3D(Type const& xCoord, Type const& yCoord, Type const& zCoord, Type const& alpha)
      : m_x(xCoord), m_y(yCoord), m_z(zCoord) {this->normalize(alpha);}


    // Copy constructor.
    template <class Type>
    Vector3D<Type>::
    Vector3D(const Vector3D<Type>& source)
      : m_x(source.m_x), m_y(source.m_y), m_z(source.m_z) {}


    // Destructor.
    template <class Type>
    Vector3D<Type>::
    ~Vector3D() {}


    // Resets the vector to (0.0, 0.0, 0.0).
    template <class Type>
    inline Vector3D<Type>&
    Vector3D<Type>::
    clear() {
      m_x = 0.0; m_y = 0.0; m_z = 0.0;
      return *this;
    }


    // Explicitly sets 3D coordinates.
    template <class Type>
    inline Vector3D<Type> const&
    Vector3D<Type>::
    setValue(Type const& xCoord, Type const& yCoord, Type const& zCoord)
    {
      m_x = xCoord; m_y = yCoord; m_z = zCoord;
      return *this;
    }


    // Explicitly sets 3D homogeneous coordinates.
    template <class Type>
    inline Vector3D<Type> const&
    Vector3D<Type>::
    setValue(Type const& xCoord, Type const& yCoord, Type const& zCoord,
             Type const& alpha)
    {
      m_x = xCoord; m_y = yCoord; m_z = zCoord; normalize(alpha);
      return *this;
    }


    // Returns the x component of the vector by value.
    template <class Type>
    inline Type const&
    Vector3D<Type>::
    x() const
    {
      return m_x;
    }


    // Returns the y component of the vector by value.
    template <class Type>
    inline Type const&
    Vector3D<Type>::
    y() const
    {
      return m_y;
    }


    // Returns the z component of the vector by value.
    template <class Type>
    inline Type const&
    Vector3D<Type>::
    z() const
    {
      return m_z;
    }


    // Assignment operator.
    template <class Type>
    Vector3D<Type>&
    Vector3D<Type>::
    operator=(const Vector3D<Type>& source)
    {
      if(&source != this) {
        setValue(source.m_x, source.m_y, source.m_z);
      }
      return *this;
    }


    // The indexing operator returns a reference to the x, y, or z
    // component of *this as if *this were a three element array.
    template <class Type>
    inline Type&
    Vector3D<Type>::
    operator[](size_t index)
    {
      switch(index) {
      case 0: return this->m_x; break;
      case 1: return this->m_y; break;
      }
      return this->m_z;
    }


    // The indexing operator returns the value of the x, y, or z
    // component of *this as if *this were a three element array.
    template <class Type>
    inline Type const&
    Vector3D<Type>::
    operator[](size_t index) const
    {
      switch(index) {
      case 0: return this->x(); break;
      case 1: return this->y(); break;
      }
      return this->z();
    }


    // Multiplies each element by a scalar.
    template <class Type>
    Vector3D<Type>&
    Vector3D<Type>::
    operator*=(Type const& scalar)
    {
      m_x *= scalar; m_y *= scalar; m_z *= scalar; return *this;
    }


    // Divides each element by a scalar.
    template <class Type>
    Vector3D<Type>&
    Vector3D<Type>::
    operator/=(Type const& scalar)
    {
      if (scalar == static_cast<Type>(0)) {
        BRICK_THROW(common::ValueException, "Vector3D::operator/=(Type)",
                    "Can't divide by zero.");
      }
      m_x /= scalar; m_y /= scalar; m_z /= scalar; return *this;
    }


    // Adds the elements of another Vector3D.
    template <class Type>
    Vector3D<Type>&
    Vector3D<Type>::
    operator+=(const Vector3D<Type>& vec)
    {
      m_x += vec.m_x; m_y += vec.m_y; m_z += vec.m_z; return *this;
    }


    // Subtracts the elements of another Vector3D.
    template <class Type>
    Vector3D<Type>&
    Vector3D<Type>::
    operator-=(const Vector3D<Type>& vec)
    {
      m_x -= vec.m_x; m_y -= vec.m_y; m_z -= vec.m_z; return *this;
    }


    // Returns a Vector3D equal to *this, but with each element negated.
    template <class Type>
    Vector3D<Type>
    Vector3D<Type>::
    operator-()
    {
      return Vector3D<Type>(-m_x, -m_y, -m_z);
    }


    // Private member functions.
    template <class Type>
    inline void
    Vector3D<Type>::
    normalize(Type const& alpha)
    {
      if(alpha == static_cast<Type>(1.0)) {return;}
      if(alpha == static_cast<Type>(0.0)) {
        BRICK_THROW(common::ValueException, "Vector3D::normalize()",
                    "Bad alpha (0.0).");
      }
      m_x /= alpha; m_y /= alpha; m_z /= alpha;
      return;
    }


    // This operator returns the elementwise sum of two Vector3D instances.
    template <class Type>
    Vector3D<Type>
    operator+(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1)
    {
      return Vector3D<Type>(vector0.x() + vector1.x(),
                      vector0.y() + vector1.y(),
                      vector0.z() + vector1.z());
    }


    // This operator returns the elementwise difference of two Vector3D
    // instances.
    template <class Type>
    Vector3D<Type>
    operator-(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1)
    {
      return Vector3D<Type>(vector0.x() - vector1.x(),
                      vector0.y() - vector1.y(),
                      vector0.z() - vector1.z());
    }


    // This operator returns the elementwise product of two Vector3D instances.
    template <class Type>
    Vector3D<Type>
    operator*(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1)
    {
      return Vector3D<Type>(vector0.x() * vector1.x(),
                      vector0.y() * vector1.y(),
                      vector0.z() * vector1.z());
    }


    // This operator returns the elementwise dividend of two Vector3D instances.
    template <class Type>
    Vector3D<Type>
    operator/(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1)
    {
      return Vector3D<Type>(vector0.x() / vector1.x(),
                      vector0.y() / vector1.y(),
                      vector0.z() / vector1.z());
    }


    // This operator adds a scalar and a Vector3D.
    template <class Type>
    Vector3D<Type> operator+(const Vector3D<Type>& vector0, Type const& scalar0)
    {
      return Vector3D<Type>(vector0.x() + scalar0,
                            vector0.y() + scalar0,
                            vector0.z() + scalar0);
    }


    // This operator subtracts a scalar from a Vector3D.
    template <class Type>
    Vector3D<Type> operator-(const Vector3D<Type>& vector0, Type const& scalar0)
    {
      return Vector3D<Type>(vector0.x() - scalar0,
                            vector0.y() - scalar0,
                            vector0.z() - scalar0);
    }


    // This operator multiplies a Vector3D by scalar.
    template <class Type>
    Vector3D<Type> operator*(const Vector3D<Type>& vector0, Type const& scalar0)
    {
      return Vector3D<Type>(vector0.x() * scalar0,
                            vector0.y() * scalar0,
                            vector0.z() * scalar0);
    }


    // This operator divides a Vector3D by scalar.
    template <class Type>
    Vector3D<Type> operator/(const Vector3D<Type>& vector0, Type const& scalar0)
    {
      return Vector3D<Type>(vector0.x() / scalar0,
                            vector0.y() / scalar0,
                            vector0.z() / scalar0);
    }


    // This operator checks the supplied vectors for equality.
    template <class Type>
    bool operator==(Vector3D<Type> const& vector0,
                    Vector3D<Type> const& vector1)
    {
      return((vector0.x() == vector1.x()) &&
             (vector0.y() == vector1.y()) &&
             (vector0.z() == vector1.z()));
    }


    // This operator checks the supplied vectors for inequality.
    template <class Type>
    bool operator!=(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1)
    {
      return(!operator==(vector0, vector1));
    }


    // This operator adds a scalar value to each element of a Vector3D
    // instance.
    template <class Type>
    Vector3D<Type> operator+(Type const& scalar0, const Vector3D<Type>& vector0)
    {
      return vector0 + scalar0;
    }


    // This operator multiplies a scalar value with each element of a
    // Vector3D instance.
    template <class Type>
    Vector3D<Type> operator*(Type const& scalar0, const Vector3D<Type>& vector0)
    {
      return vector0 * scalar0;
    }


    // This function outputs a text representation of a Vector3D
    // instance to a std::ostream.
    template <class Type>
    std::ostream& operator<<(std::ostream& stream, const Vector3D<Type>& vector0)
    {
      stream << "Vector3D(" << vector0.x() << ", " << vector0.y() << ", "
             << vector0.z() << ")";
      return stream;
    }


    // This function sets the value of a Vector3D instance from a
    // std::istream.
    template <class Type>
    std::istream& operator>>(std::istream& stream, Vector3D<Type>& vector0)
    {
      const char intro[] = "Vector3D(";
      const char intermission[] = ",";
      const char outro[] = ")";
      Type x, y, z;
      char inChar;
      size_t index;

      for(index = 0; index < strlen(intro); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != intro[index]) {
          // Hmm.  g++ diverges from the standard?
          // stream.clear(ios_base::failbit);
          stream.clear(std::ios::failbit);
          return stream;
        }
      }
      stream >> x;
      for(index = 0; index < strlen(intermission); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != intermission[index]) {
          // stream.clear(ios_base::failbit);
          stream.clear(std::ios::failbit);
          return stream;
        }
      }
      stream >> y;
      for(index = 0; index < strlen(intermission); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != intermission[index]) {
          // stream.clear(ios_base::failbit);
          stream.clear(std::ios::failbit);
          return stream;
        }
      }
      stream >> z;
      for(index = 0; index < strlen(outro); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != outro[index]) {
          // stream.clear(ios_base::failbit);
          stream.clear(std::ios::failbit);
          return stream;
        }
      }
      if(stream) {
        vector0.setValue(x, y, z);
      }
      return stream;
    }

  } // namespace numeric

} // namespace brick

#endif  /* #ifndef BRICK_NUMERIC_VECTOR3D_IMPL_HH */
