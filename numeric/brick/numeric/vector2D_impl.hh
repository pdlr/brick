/**
***************************************************************************
* @file brick/numeric/vector2D_impl.hh
*
* Implementation file for Vector2D class template.
*
* Copyright (C) 2000-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_VECTOR2D_IMPL_HH
#define BRICK_NUMERIC_VECTOR2D_IMPL_HH

// This file is included by vector2D.hh, and should not be directly included
// by user code, so no need to include vector2D.hh here.
// 
// #include <brick/numeric/vector2D.hh>

namespace brick {

  namespace numeric {

    // The default constructor initializes to (0, 0).
    template <class Type>
    Vector2D<Type>::
    Vector2D()
      : m_x(0),
        m_y(0)
    {
      // Empty.
    }
    

    // This constructor explicitly sets the 2D coordinates.
    template <class Type>
    Vector2D<Type>::
    Vector2D(Type const& xCoord, Type const& yCoord)
      : m_x(xCoord),
        m_y(yCoord)
    {
      // Empty.
    }
        


    // This constructor explicitly sets 2D homogeneous coordinates.
    template <class Type>
    Vector2D<Type>::
    Vector2D(Type const& xCoord, Type const& yCoord, Type const& alpha)
      : m_x(xCoord),
        m_y(yCoord)
    {
      this->normalize(alpha);
    }


    // The copy constructor deep copies its argument.
    template <class Type>
    Vector2D<Type>::
    Vector2D(const Vector2D<Type>& vec)
      : m_x(vec.m_x),
        m_y(vec.m_y)
    {
      // Empty.
    }


    // The destructor destroys the Vector2D instance.
    template <class Type>
    Vector2D<Type>::
    ~Vector2D()
    {
      // Empty.
    }


    // Resets the vector to (0.0, 0.0).
    template <class Type>
    inline Vector2D<Type>&
    Vector2D<Type>::
    clear()
    {
      m_x = 0.0; m_y = 0.0;
      return *this;
    }

      
    // This member function returns the X component of the Vector2D.
    template <class Type>
    inline Type const&
    Vector2D<Type>::
    getX() const
    {
      return m_x;
    }


    // This member function returns the Y component of the Vector2D.
    template <class Type>
    inline Type const&
    Vector2D<Type>::
    getY() const
    {
      return m_y;
    }

    
    // This member function sets the X component of the Vector2D.
    template <class Type>
    inline Type const&
    Vector2D<Type>::
    setX(Type const& newX)
    {
      m_x = newX;
      return m_x;
    }


    // This member function sets the Y component of the Vector2D.
    template <class Type>
    inline Type const&
    Vector2D<Type>::
    setY(Type const& newY)
    {
      m_y = newY;
      return m_y;
    }


    // This member function explicitly sets the sets coordinates of
    // the Vector2D instance.
    template <class Type>
    inline Vector2D<Type> const&
    Vector2D<Type>::
    setValue(Type const& xCoord, Type const& yCoord)
    {
      m_x = xCoord;
      m_y = yCoord;
      return *this;
    }


    // This member function explicitly sets 2D homogeneous
    // coordinates.
    template <class Type>
    inline Vector2D<Type> const&
    Vector2D<Type>::
    setValue(Type const& xCoord, Type const& yCoord, Type const& alpha)
    {
      m_x = xCoord;
      m_y = yCoord;
      this->normalize(alpha);
      return *this;
    }


#if 0
    
    // This member function returns the X component of the Vector2D
    // by value.
    template <class Type>
    inline Type
    Vector2D<Type>::
    x() const
    {
      return m_x;
    }


    // This member function returns the Y component of the Vector2D
    // by value.
    template <class Type>
    inline double
    Vector2D<Type>::
    y() const
    {
      return m_y;
    }
    
#endif /* #if 0 */

    // The assignment operator deep copies its argument.
    template <class Type>
    Vector2D<Type>&
    Vector2D<Type>::
    operator=(const Vector2D<Type>& vec)
    {
      // Self-assignment OK.
      setValue(vec.m_x, vec.m_y);
      return *this;
    }


    // The indexing operator returns a reference to the x, or y
    // component of *this as if *this were a two element array.
    template <class Type>
    inline Type&
    Vector2D<Type>::
    operator[](size_t index)
    {
      if(index == 0) {
        return this->m_x;
      }
      return this->m_y;
    }

          
    // The indexing operator returns the value of the x, or y
    // component of *this as if *this were a two element array.
    template <class Type>
    inline Type
    Vector2D<Type>::
    operator[](size_t index) const
    {
      if(index == 0) {return this->x();}
      return this->y();
    }
          

    // This operator multiplies each component of the Vector2D instance
    // by a scalar.
    template <class Type>
    Vector2D<Type>&
    Vector2D<Type>::
    operator*=(Type const& scalar)
    {
      m_x *= scalar; m_y *= scalar; return *this;
    }


    // This operator divides each component of the Vector2D instance
    // by a scalar.
    template <class Type>
    Vector2D<Type>&
    Vector2D<Type>::
    operator/=(Type const& scalar)
    {
      if (scalar == 0.0) {
        BRICK_THROW(common::ValueException, "Vector2D::operator/=()",
                    "Bad scalar: Divide by Zero\n");
      }
      m_x /= scalar; m_y /= scalar; return *this;
    }


    // This operator adds a scalar to each component of the Vector2D
    // instance.
    template <class Type>
    Vector2D<Type>&
    Vector2D<Type>::
    operator+=(const Vector2D<Type>& vec)
    {
      m_x += vec.m_x; m_y += vec.m_y; return *this;
    }


    // This operator subtracts a scalar from each component of the
    // Vector2D instance.
    template <class Type>
    Vector2D<Type>&
    Vector2D<Type>::
    operator-=(const Vector2D<Type>& vec)
    {
      m_x -= vec.m_x; m_y -= vec.m_y; return *this;
    }


    // This operator returns a Vector2D equal to *this, but with each
    // element negated.
    template <class Type>
    Vector2D<Type>
    Vector2D<Type>::
    operator-()
    {
      return Vector2D<Type>(-m_x, -m_y);
    }


    template <class Type>
    inline void
    Vector2D<Type>::
    normalize(Type const& alpha)
    {
      if(alpha == static_cast<Type>(1.0)) {return;}
      if(alpha == static_cast<Type>(0.0)) {
        BRICK_THROW(common::ValueException, "Vector2D::normalize()",
                    "Bad alpha (0.0).");
      }
      m_x /= alpha; m_y /= alpha;
      return;
    }

    template <class Type>
    inline Vector2D<Type>
    operator+(Type const& scalar0, const Vector2D<Type>& vector0)
    {
      return vector0 + scalar0;
    }
  
    template <class Type>
    inline Vector2D<Type>
    operator*(Type const& scalar0, const Vector2D<Type>& vector0)
    {
      return vector0 * scalar0;
    }

      
    template <class Type>
    Vector2D<Type>
    operator+(const Vector2D<Type>& vector0, const Vector2D<Type>& vector1)
    {
      return Vector2D<Type>(vector0.x() + vector1.x(),
                            vector0.y() + vector1.y());
    }
  
    template <class Type>
    Vector2D<Type>
    operator-(const Vector2D<Type>& vector0, const Vector2D<Type>& vector1)
    {
      return Vector2D<Type>(vector0.x() - vector1.x(),
                            vector0.y() - vector1.y());
    }
  
    template <class Type>
    Vector2D<Type>
    operator*(const Vector2D<Type>& vector0, const Vector2D<Type>& vector1)
    {
      return Vector2D<Type>(vector0.x() * vector1.x(),
                            vector0.y() * vector1.y());
    }
  
    template <class Type>
    Vector2D<Type>
    operator/(const Vector2D<Type>& vector0, const Vector2D<Type>& vector1)
    {
      return Vector2D<Type>(vector0.x() / vector1.x(),
                            vector0.y() / vector1.y());
    }
  
    template <class Type>
    Vector2D<Type>
    operator+(const Vector2D<Type>& vector0, Type const& scalar)
    {
      return Vector2D<Type>(vector0.x() + scalar,
                            vector0.y() + scalar);
    }

    template <class Type>
    Vector2D<Type>
    operator-(const Vector2D<Type>& vector0, Type const& scalar)
    {
      return Vector2D<Type>(vector0.x() - scalar,
                            vector0.y() - scalar);
    }

    template <class Type>
    Vector2D<Type>
    operator*(const Vector2D<Type>& vector0, Type const& scalar)
    {
      return Vector2D<Type>(vector0.x() * scalar,
                            vector0.y() * scalar);
    }

    template <class Type>
    Vector2D<Type>
    operator/(const Vector2D<Type>& vector0, Type const& scalar)
    {
      return Vector2D<Type>(vector0.x() / scalar,
                            vector0.y() / scalar);
    }

    template <class Type>
    bool
    operator==(const Vector2D<Type>& vector0, const Vector2D<Type>& vector1)
    {
      return((vector0.x() == vector1.x()) &&
             (vector0.y() == vector1.y()));
    }

    template <class Type>
    bool
    operator!=(const Vector2D<Type>& vector0, const Vector2D<Type>& vector1)
    {
      return(!operator==(vector0, vector1));
    }

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Vector2D<Type>& vector0)
    {
      stream << "Vector2D(" << vector0.x() << ", " << vector0.y() << ")";
      return stream;
    }

    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Vector2D<Type>& vector0)
    {
      const char intro[] = "Vector2D<Type>(";
      const char intermission[] = ",";
      const char outro[] = ")";
      Type x, y;
      char inChar;
      size_t index;

      for(index = 0; index < strlen(intro); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != intro[index]) {
          stream.clear(std::ios_base::failbit);
          return stream;
        }
      }
      stream >> x;
      for(index = 0; index < strlen(intermission); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != intermission[index]) {
          stream.clear(std::ios_base::failbit);
          return stream;
        }
      }
      stream >> y;
      for(index = 0; index < strlen(outro); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != outro[index]) {
          stream.clear(std::ios_base::failbit);
          return stream;
        }
      }
      if(stream) {
        vector0.setValue(x, y);
      }
      return stream;
    }

  } // namespace numeric

} // namespace brick

#endif  /* #ifndef BRICK_NUMERIC_VECTOR2D_IMPL_HH */
