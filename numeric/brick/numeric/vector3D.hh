/**
***************************************************************************
* @file brick/numeric/vector3D.hh
*
* Header file declaring Vector3D class template.
*
* Copyright (C) 2000-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_VECTOR3D_HH
#define BRICK_NUMERIC_VECTOR3D_HH

#include <iostream>
#include <brick/common/exception.hh>

namespace brick {

  namespace numeric {
    
    /**
     ** The Vector3D class template represents a 3D vector in which
     ** the elements are of a user specified type.
     **/
    template <class Type>
    class Vector3D {
    public:
      /** 
       *  Default constructor initializes to (0, 0, 0).
       */
      inline
      Vector3D();

      /** 
       * Explicitly sets 3D coordinates.
       * 
       * @param xCoord After construction, the vector will have this value as
       * its X coordinate.
       * @param yCoord After construction, the vector will have this value as
       * its Y coordinate.
       * @param zCoord After construction, the vector will have this value as
       * its Z coordinate.
       */
      inline
      Vector3D(Type const& xCoord, Type const& yCoord, Type const& zCoord);

      
      /** 
       * Explicitly sets 3D homogeneous coordinates.  A 3D homogeneous
       * vector has the form (xCoord, yCoord, zCoord, alpha), and
       * corresponds to the 3D point (xCoord/alpha, yCoord/alpha,
       * zCoord/alpha).
       * 
       * @param xCoord The homogeneous X coordinate.
       * @param yCoord The homogeneous Y coordinate.
       * @param zCoord The homogeneous Z coordinate.
       * @param alpha Scale factor.
       */
      inline
      Vector3D(Type const& xCoord, Type const& yCoord, Type const& zCoord,
               Type const& alpha);

      
      /** 
       * Copy constructor.
       * 
       * @param source The Vector3D to be copied.
       */
      inline
      Vector3D(const Vector3D<Type>& source);


      /** 
       * Destructor.
       */
      inline
      ~Vector3D();

      
      /** 
       * Resets the vector to (0.0, 0.0, 0.0).
       *
       * @return A reference to *this.
       */
      inline Vector3D<Type>&
      clear();

      
      /** 
       * This member function returns the X component of the Vector2D.
       * 
       * @return The return value is the X coordinate.
       */
      inline Type const&
      getX() const {return m_x;}


      /** 
       * This member function returns the Y component of the Vector2D.
       * 
       * @return The return value is the Y coordinate.
       */
      inline Type const&
      getY() const {return m_y;}

    
      /** 
       * This member function returns the Y component of the Vector2D.
       * 
       * @return The return value is the Y coordinate.
       */
      inline Type const&
      getZ() const {return m_z;}

    
      /** 
       * Explicitly sets 3D coordinates.
       * 
       * @param xCoord The desired X coordinate.
       * @param yCoord The desired Y coordinate.
       * @param zCoord The desired Z coordinate.
       */
      inline void
      setValue(Type const& xCoord, Type const& yCoord, Type const& zCoord);


      /** 
       * This member function sets the X component of the Vector2D.
       *
       * @param newX This parameter is the value to which the X
       * coordinate will be set.
       *
       * @return The return value is the new X coordinate.
       */
      inline Type const&
      setX(Type const& newX) {m_x = newX; return m_x;}


      /** 
       * This member function sets the Y component of the Vector2D.
       *
       * @param newY This parameter is the value to which the Y
       * coordinate will be set.
       *
       * @return The return value is the new Y coordinate.
       */
      inline Type const&
      setY(Type const& newY) {m_y = newY; return m_y;}


      /** 
       * This member function sets the Y component of the Vector2D.
       *
       * @param newY This parameter is the value to which the Y
       * coordinate will be set.
       *
       * @return The return value is the new Y coordinate.
       */
      inline Type const&
      setZ(Type const& newZ) {m_z = newZ; return m_z;}


      /** 
       * Explicitly sets 3D homogeneous coordinates.
       * 
       * @param xCoord The homogeneous X coordinate.
       * @param yCoord The homogeneous Y coordinate.
       * @param zCoord The homogeneous Z coordinate.
       * @param alpha Scale factor.
       */
      inline void setValue(Type const& xCoord,
                           Type const& yCoord,
                           Type const& zCoord,
                           Type const& alpha);

      
#if 0
      /** 
       * Returns the x component of the vector by reference.
       * 
       * @return A reference to the x component of the vector.
       */
      inline Type& x() {return m_x;}
#endif /* #if 0 */


      /** 
       * Returns the x component of the vector by value.
       * 
       * @return The value of the x component of the vector.
       */
      inline Type const& x() const;

#if 0
      /** 
       * Returns the y component of the vector by reference.
       * 
       * @return A reference to the y component of the vector.
       */
      inline Type& y() {return m_y;}
#endif /* #if 0 */

      /** 
       * Returns the y component of the vector by value.
       * 
       * @return The value of the y component of the vector.
       */
      inline Type const& y() const;
    
#if 0
      /** 
       * Returns the z component of the vector by reference.
       * 
       * @return A reference to the z component of the vector.
       */
      inline Type& z() {return m_z;}
#endif /* #if 0 */

      
      /** 
       * Returns the z component of the vector by value.
       * 
       * @return The value of the z component of the vector.
       */
      inline Type const& z() const;

      
      /**
       * Assignment operator.
       * 
       * @param source The vector to be copied.
       * @return Reference to *this.
       */
      inline Vector3D<Type>&
      operator=(const Vector3D<Type>& source);
      

      /** 
       * The indexing operator returns a reference to the x, y, or z
       * component of *this as if *this were a three element array.
       * Out of bounds indices will return this->z().
       * 
       * @param index This argument is the index into *this.
       * 
       * @return The return value is the selected component of *this.
       */
      inline Type&
      operator[](size_t index);

          
      /** 
       * The indexing operator returns the value of the x, y, or z
       * component of *this as if *this were a three element array.
       * Out of bounds indices will return this-z().
       * 
       * @param index This argument is the index into *this.
       * 
       * @return The return value is the selected component of *this.
       */
      inline Type const& operator[](size_t index) const;

          
      /** 
       * Multiplies each element by a scalar.
       * 
       * @param scalar X, Y, and Z values will be multiplied by this value.
       * @return Reference to *this.
       */
      inline Vector3D<Type>&
      operator*=(Type const& scalar);

      
      /** 
       * Divides each element by a scalar.
       * 
       * @param scalar X, Y, and Z values will be divided by this value.
       * @return Reference to *this.
       */
      inline Vector3D<Type>&
      operator/=(Type const& scalar);


      /** 
       * Adds the elements of another Vector3D.
       * 
       * @param vec The elements of vec will be added to *this.
       * @return Reference to *this.
       */
      inline Vector3D<Type>&
      operator+=(const Vector3D<Type>& vec);


      /** 
       * Subtracts the elements of another Vector3D.
       * 
       * @param vec The elements of vec will be subtracted from *this.
       * @return Reference to *this.
       */
      inline Vector3D<Type>&
      operator-=(const Vector3D<Type>& vec);


      /** 
       * Returns a Vector3D equal to *this, but with each element negated.
       * 
       * @return The result of the negation.
       */
      inline Vector3D<Type>
      operator-();

    private:
      // Private member functions.
      void
      normalize(Type const& alpha);

      // Private data members.
      Type m_x;
      Type m_y;
      Type m_z;
    }; // class Vector3D


    /* ============== Non-member function declarations ============== */
  
    /** 
     * This operator returns the elementwise sum of two Vector3D instances.
     * 
     * @param vector0 This is the first of the two Vector3D instances to
     * be added.
     * @param vector1 This is the second of the two Vector3D instances to
     * be added.
     * @return A Vector3D instance in which the value of each element is
     * equal to the sum of the corresponding elements of the two arguments.
     */
    template <class Type>
    Vector3D<Type>
    operator+(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1);

    
    /** 
     * This operator returns the elementwise difference of two Vector3D
     * instances.
     * 
     * @param vector0 This is the first of the two Vector3D instances to
     * be subtracted.
     * @param vector1 This is the second of the two Vector3D instances to
     * be subtracted.
     * @return A Vector3D instance in which the value of each element is
     * equal to the difference of the corresponding elements of the two
     * arguments.
     */
    template <class Type>
    Vector3D<Type>
    operator-(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1);

    
    /** 
     * This operator returns the elementwise product of two Vector3D instances.
     * 
     * @param vector0 This is the first of the two Vector3D instances to
     * be multiplied.
     * @param vector1 This is the second of the two Vector3D instances to
     * be multiplied.
     * @return A Vector3D instance in which the value of each element is
     * equal to the product of the corresponding elements of the two arguments.
     */
    template <class Type>
    Vector3D<Type>
    operator*(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1);

    
    /** 
     * This operator returns the elementwise dividend of two Vector3D instances.
     * 
     * @param vector0 This is the Vector3D instance whose element values
     * are to be divided.
     * @param vector1 This is the Vector3D instance by whose elements
     * the first argument is to be divided.
     * @return A Vector3D instance in which the value of each element is
     * equal to the corresponding value of the first argument divided by
     * the corresponding value of the second argument.
     */
    template <class Type>
    Vector3D<Type>
    operator/(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1);

    
    /** 
     * This operator adds a scalar and a Vector3D.
     * 
     * @param vector0 This is the Vector3D instance to which the scalar
     * should be added.
     * @param scalar0 This is amount which should be added to each
     * element of argument vector0.
     * @return A Vector3D instance in which the value of each element is
     * equal to the corresponding value of the first argument plus the
     * value of the second argument.
     */
    template <class Type>
    Vector3D<Type>
    operator+(const Vector3D<Type>& vector0, Type const& scalar0);

    
    /** 
     * This operator subtracts a scalar from a Vector3D.
     * 
     * @param vector0 This is the Vector3D instance from which the scalar
     * should be subtracted.
     * @param scalar0 This is amount which should be subtracted from each
     * element of argument vector0.
     * @return A Vector3D instance in which the value of each element is
     * equal to the corresponding value of the first argument minus the
     * value of the second argument.
     */
    template <class Type>
    Vector3D<Type>
    operator-(const Vector3D<Type>& vector0, Type const& scalar0);

    
    /** 
     * This operator multiplies a Vector3D by scalar.
     * 
     * @param vector0 This is the Vector3D instance which is to be
     * multiplied by the scalar.
     * @param scalar0 This is amount by which should argument vector0 is
     * to be multiplied.
     * @return A Vector3D instance in which the value of each element is
     * equal to the corresponding value of the first argument multiplied by
     * the value of the second argument.
     */
    template <class Type>
    Vector3D<Type>
    operator*(const Vector3D<Type>& vector0, Type const& scalar0);

    
    /** 
     * This operator divides a Vector3D by scalar.
     * 
     * @param vector0 This is the Vector3D instance which is to be
     * divided by the scalar.
     * @param scalar0 This is amount by which should argument vector0 is
     * to be divided.
     * @return A Vector3D instance in which the value of each element is
     * equal to the corresponding value of the first argument divided by
     * the value of the second argument.
     */
    template <class Type>
    Vector3D<Type>
    operator/(const Vector3D<Type>& vector0, Type const& scalar0);

    
    /**
     * This operator checks the supplied vectors for equality.
     *
     * @param  vector0  First vector to compare.
     * @param  vector1  Second vector to compare.
     * @return  Result of comparing @p vector0 to @p vector1 for equality.
     */
    template <class Type>
    bool
    operator==(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1);

    
    /**
     * This operator checks the supplied vectors for inequality.
     *
     * @param  vector0  First vector to compare.
     * @param  vector1  Second vector to compare.
     * @return  Result of comparing @p vector0 to @p vector1 for
     *          inequality.
     */
    template <class Type>
    bool
    operator!=(const Vector3D<Type>& vector0, const Vector3D<Type>& vector1);


    /** 
     * This operator adds a scalar value to each element of a Vector3D
     * instance.
     * 
     * @param scalar0 Scalar argument of the addition.
     *
     * @param vector0 Vector argument of the addition.
     *
     * @return Vector3D instance in which the value of each element is
     * the sum of the scalar argument and the corresponding element of
     * the Vector3D argument.
     */
    template <class Type>
    Vector3D<Type>
    operator+(Type const& scalar0, const Vector3D<Type>& vector0);


    /** 
     * This operator multiplies a scalar value with each element of a
     * Vector3D instance.
     * 
     * @param scalar0 Scalar argument of the multiplication.
     *
     * @param vector0 Vector argument of the multiplication.
     *
     * @return Vector3D instance in which the value of each element is
     * the product of the scalar argument and the corresponding element
     * of the Vector3D argument.
     */
    template <class Type>
    Vector3D<Type>
    operator*(Type const& scalar0, const Vector3D<Type>& vector0);


    /** 
     * This function outputs a text representation of a Vector3D
     * instance to a std::ostream.  The output format looks like this:
     *
     * Vector3D(1.0, 2.0, 3.0)
     *
     * @param stream This argument is a reference to the the output
     * stream.
     *
     * @param vector0 This argument is a const reference to the
     * Vector3D instance to be output.
     *
     * @return The return value is a reference to the input stream after
     * the write has taken place.
     */
    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Vector3D<Type>& vector0);


    /** 
     * This function sets the value of a Vector3D instance from a
     * std::istream.  The input format is as described for
     * operator<<(std::ostream&, const Vector3D&) above.
     * 
     * @param stream This argument is a reference to the the input
     * stream from which to read.
     *
     * @param vector0 This argument is a reference to the Vector3D
     * which will take the input.
     *
     * @return The return value is a reference to the input stream after
     * the read has taken place.
     */
    template <class Type>
    std::istream&
    operator>>(std::istream& stream, Vector3D<Type>& vector0);

  } // namespace numeric

} // namespace brick


// Definitions of inline and template functions.
#include <brick/numeric/vector3D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_VECTOR3D_HH */
