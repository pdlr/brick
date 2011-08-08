/**
***************************************************************************
* @file index3D.cpp
*
* Source file defining Index3D class.
*
* Copyright (C) 2000-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/numeric/index3D.hh>

namespace brick {

  namespace numeric {
    
    Index3D operator+(const Index3D& index0, const Index3D& index1)
    {
      return Index3D(index0.getRow() + index1.getRow(),
                     index0.getColumn() + index1.getColumn(),
                     index0.getSlice() + index1.getSlice());
    }

  
    Index3D operator-(const Index3D& index0, const Index3D& index1)
    {
      return Index3D(index0.getRow() - index1.getRow(),
                     index0.getColumn() - index1.getColumn(),
                     index0.getSlice() - index1.getSlice());
    }

  
    Index3D operator*(const Index3D& index0, const Index3D& index1)
    {
      return Index3D(index0.getRow() * index1.getRow(),
                     index0.getColumn() * index1.getColumn(),
                     index0.getSlice() * index1.getSlice());
    }

  
    Index3D operator/(const Index3D& index0, const Index3D& index1)
    {
      return Index3D(index0.getRow() / index1.getRow(),
                     index0.getColumn() / index1.getColumn(),
                     index0.getSlice() / index1.getSlice());
    }

  
    Index3D operator+(const Index3D& index0, int scalar)
    {
      return Index3D(index0.getRow() + scalar,
                     index0.getColumn() + scalar,
                     index0.getSlice() + scalar);
    }

  
    Index3D operator-(const Index3D& index0, int scalar)
    {
      return Index3D(index0.getRow() - scalar,
                     index0.getColumn() - scalar,
                     index0.getSlice() - scalar);
    }

  
    Index3D operator*(const Index3D& index0, int scalar)
    {
      return Index3D(index0.getRow() * scalar,
                     index0.getColumn() * scalar,
                     index0.getSlice() * scalar);
    }

  
    Index3D operator/(const Index3D& index0, int scalar)
    {
      return Index3D(index0.getRow() / scalar,
                     index0.getColumn() / scalar,
                     index0.getSlice() / scalar);
    }

  
    bool
    operator==(const Index3D& index0, const Index3D& index1)
    {
      return((index0.getRow() == index1.getRow()) &&
             (index0.getColumn() == index1.getColumn()) &&
             (index0.getSlice() == index1.getSlice()));
    }

  
    bool
    operator!=(const Index3D& index0, const Index3D& index1)
    {
      return(!operator==(index0, index1));
    }

  
    std::ostream&
    operator<<(std::ostream& stream, const Index3D& index0)
    {
      stream << "Index3D(" << index0.getSlice() << ", "
             << index0.getRow() << ", "
             << index0.getColumn() << ")";
      return stream;
    }

  
    std::istream&
    operator>>(std::istream& stream, Index3D& index0)
    {
      const char intro[] = "Index3D(";
      const char intermission[] = ",";
      const char outro[] = ")";
      int row, column, slice;
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
      stream >> slice;
      for(index = 0; index < strlen(intermission); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != intermission[index]) {
          stream.clear(std::ios_base::failbit);
          return stream;
        }
      }
      stream >> row;
      for(index = 0; index < strlen(intermission); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != intermission[index]) {
          stream.clear(std::ios_base::failbit);
          return stream;
        }
      }
      stream >> column;
      for(index = 0; index < strlen(outro); ++index) {
        inChar = 0;
        stream >> inChar;
        if(inChar != outro[index]) {
          stream.clear(std::ios_base::failbit);
          return stream;
        }
      }
      if(stream) {
        index0.setValue(slice, row, column);
      }
      return stream;
    }

  } // namespace numeric

} // namespace brick
