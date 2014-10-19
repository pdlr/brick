/**
***************************************************************************
* @file brick/common/versionNumber.hh>
*
* Header file declaring VersionNumber class.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_VERSIONNUMBER_HH
#define BRICK_COMMON_VERSIONNUMBER_HH

namespace brick {

  namespace common {
    
    /**
     ** The VersionNumber class permits convenient input, output, and
     ** comparison of version numbers with formats like "1.23" ,
     ** "3.2.55" , "6.5.4.3.2" , etc.
     **/
    class VersionNumber
    {
    public:

      /** 
       * The default constructor initializes to version 0.
       */
      VersionNumber(UInt32 numberOfFieldsHint = 1) 
	: m_fields(numberOfFieldsHint, 0) {}

      
      /** 
       * This constructor creates a VersionNumber from a string, like "1.2.3".
       * 
       * @param version A string representing the version number.  It
       * must be a series of integers seperated by "." characters.
       */
      VersionNumber(std::string const& versionString, 
		    UInt32 numberOfFieldsHint = 1) 
	: m_fields(numberOfFieldsHint) {
	std::istringstream buffer(versionString);
	buffer >> *this;
      }

			  
      /** 
       * Copy constructor.
       * 
       * @param source The VersionNumber instance to be copied.
       */
      VersionNumber(const VersionNumber<Type>& source)
        : m_fields(source.m_fields) {}

      
      /** 
       * Destroys the VersionNumber instance.
       */
      ~VersionNumber() {}

      
      /** 
       * The pre-increment operator actually adds the value of stride
       * (which was set at construction time) to *this.  That is, if the
       * value of stride is N, then incrementing a VersionNumber
       * instance makes it point to an array element N spaces away.
       * 
       * @return A Reference to the incremented VersionNumber instance.
       */
      VersionNumber<Type>& operator>>() {
	m_fields.clear();
	std::vector<Int32>::back_inserter iter(m_fields);
	buffer >> (*iter);
	++iter;
	while(iter != m_fields.end()) {
	  buffer >> Expect('.') >> iter;
	  ++iter;
	}
	if(!buffer) {
	  std::ostringstream message;
	  message << "Invalid input string: \"" << versionString << "\"";
	  BRICK_THROW(ValueException, "VersionNumber::VersionNumber()",
		      message.str().c_str());
	}
      }

m_pointer += m_stride; return *this;}

      
      /** 
       * The post-increment operator is just like the pre-increment
       * operator, above, but returns a copy of the VersionNumber which
       * has not been incremented.
       * 
       * @return An un-incremented copy of *this.
       */
      VersionNumber<Type> operator++(int) {
        VersionNumber<Type> result(*this);
        ++(*this);
        return result;
      }

      
      /** 
       * The pre-decrement operator is just like the pre-increment
       * operator, above, but decrements instead of incrementing.
       * 
       * @return A Reference to the decremented VersionNumber instance.
       */
      VersionNumber<Type>& operator--() {m_pointer -= m_stride; return *this;}


      /** 
       * The post-decrement operator is just like the post-increment
       * operator, above, but decrements instead of incrementing.
       * 
       * @return An un-decremented copy of *this.
       */
      VersionNumber<Type> operator--(int) {
        VersionNumber<Type> result(*this);
        --(*this);
        return result;
      }

      
      /** 
       * This operator has the same effect as incrementing a copy of
       * *this offset times, but runs in O(1) time.  That is, a call
       * to operator+() returns a VersionNumber instance which
       * references a place in memory that is (stride * offset) items
       * advanced from the value referenced by *this.
       * 
       * @param offset This argument specifies how many steps to
       * advance.
       * @return A copy of *this that references the new location.
       */
      VersionNumber<Type> operator+(ptrdiff_t offset) {
        VersionNumber<Type> result(*this);
        result += offset;
        return result;
      }

      
      /** 
       * This operator works just like operator+() above, but decreases
       * the VersionNumber value instead of increasing it.
       * 
       * @param offset This argument specifies how many steps to
       * move back.
       * @return A copy of *this which references the new location.
       */
      VersionNumber<Type> operator-(ptrdiff_t offset) {
        VersionNumber<Type> result(*this);
        result -= offset;
        return result;
      }

      
      /** 
       * This operator attempts to return a value N so that *(*this)
       * == *(other + N).  This value is meaningless if the two
       * VersionNumber instances have different stride values.  Note
       * that if stride is not equal to 1, it is very possible to have
       * two VersionNumber instances which are "out of phase" so that
       * there is no value of N which satisfies the above equation.
       * 
       * @param other The VersionNumber instance to subtract from
       * *this.
       * @return A difference value as described above.
       */
      ptrdiff_t operator-(const VersionNumber<Type>& other) {
        ptrdiff_t distance = m_pointer - other.m_pointer;
        return distance / m_stride;
      }

      
      /** 
       * This operator has the same effect as incrementing *this offset
       * times, but runs in O(1) time.
       * 
       * @param offset This argument specifies how many steps to
       * advance the VersionNumber.
       * @return A reference to *this.
       */
      VersionNumber<Type>& operator+=(ptrdiff_t offset) {
        m_pointer += offset * m_stride;
        return *this;
      }

      
      /** 
       * This operator has the same effect as decrementing *this offset
       * times, but runs in O(1) time.
       * 
       * @param offset This argument specifies how many steps to
       * decrement the VersionNumber.
       * @return A reference to *this.
       */
      VersionNumber<Type>& operator-=(ptrdiff_t offset) {
        m_pointer -= offset * m_stride;
        return *this;
      }

      
      /** 
       * This operator returns true if its argument references the
       * same memory location as *this, false otherwise.  Note that it
       * does not require the two VersionNumber instances to have the
       * same stride.  It is sufficient that they reference the same
       * memory location.
       * 
       * @param other A VersionNumber instance to be compared with *this.
       * @return true if other references same memory location as *this,
       * false otherwise.
       * 
       * @return The return value is true if the two VersionNumber
       * instances point to the same memory location, false otherwise.
       */
      bool operator==(const VersionNumber<Type>& other) const {
        return m_pointer == other.m_pointer;
      }

      
      /** 
       * This operator returns false if its argument references the same
       * memory location as *this, true otherwise.
       * 
       * @param other A VersionNumber instance to be compared with *this.
       * @return false if other references same memory location as *this,
       * true otherwise.
       */
      bool operator!=(const VersionNumber<Type>& other) const {
        return m_pointer != other.m_pointer;
      }

      
      /** 
       * This operator return true if the address of the memory location
       * referenced by *this is less than the address of the memory
       * location referenced by other, false otherwise.
       * 
       * @param other The VersionNumber instance to be compared with
       * *this.
       * @return true if the address of the memory location referenced
       * by *this is less than the address of the memory location
       * referenced by other, false otherwise.
       */
      bool operator<(const VersionNumber<Type>& other) const {
        return m_pointer < other.m_pointer;
      }

      
      /** 
       * This operator return true if the address of the memory location
       * referenced by *this is greater than the address of the memory
       * location referenced by other, false otherwise.
       * 
       * @param other The VersionNumber instance to be compared with
       * *this.
       * @return true if the address of the memory location referenced
       * by *this is greater than the address of the memory location
       * referenced by other, false otherwise.
       */
      bool operator>(const VersionNumber<Type>& other) const {
        return m_pointer > other.m_pointer;
      }

      
      /** 
       * This operator return true if the address of the memory location
       * referenced by *this is less than or equal to the address of the
       * memory location referenced by other, false otherwise.
       * 
       * @param other The VersionNumber instance to be compared with
       * *this.
       * @return true if the address of the memory location referenced
       * by *this is less than or equal to the address of the memory
       * location referenced by other, false otherwise.
       */
      bool operator<=(const VersionNumber<Type>& other) const {
        return m_pointer <= other.m_pointer;
      }

      
      /** 
       * This operator return true if the address of the memory location
       * referenced by *this is greater than or equal to the address of
       * the memory location referenced by other, false otherwise.
       * 
       * @param other The VersionNumber instance to be compared with
       * *this.
       * @return true if the address of the memory location referenced
       * by *this is greater than or equal to the address of the memory
       * location referenced by other, false otherwise.
       */
      bool operator>=(const VersionNumber<Type>& other) const {
        return m_pointer >= other.m_pointer;
      }

      
      /** 
       * The assignment operator copies both address and stride from its
       * argument.
       * 
       * @param source The VersionNumber instance to be copied.
       * @return A reference to *this.
       */
      VersionNumber<Type>& operator=(const VersionNumber<Type>& source) {
        m_pointer = source.m_pointer;
        m_stride = source.m_stride;
        return *this;
      }

    private:
      Type* m_pointer;
      int m_stride;
    };

  } // namespace common
  
} // namespace brick

#endif /* #ifndef BRICK_COMMON_VERSIONNUMBER_HH */
