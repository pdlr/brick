/**
***************************************************************************
* @file brick/computerVision/kernel.hh
*
* Header file declaring Kernel class.
*
* Copyright (C) 2006, 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KERNEL_HH
#define BRICK_COMPUTERVISION_KERNEL_HH

#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This class template represents a 2D convolution kernel.  The template
     ** parameter indicates the type of the element in the kernel.
     **/
    template <class ELEMENT_TYPE>
    class Kernel
    {
    public:

      /* ******** Public typedefs ******** */

      typedef ELEMENT_TYPE ElementType;


      /* ******** Public member functions ******** */

      /**
       * Default constructor initializes to zero size.
       */
      Kernel();


      /**
       * The copy constructor does a deep copy.
       *
       * @param source The Kernel instance to be copied.
       */
      Kernel(const Kernel<ElementType>& source);


      /**
       * This constructor allows us to implicitly make a Kernel
       * instance from an Array2D.  As with the copy constructor, the
       * newly created kernel has its own copy of the data in the copied
       * array.
       *
       * @param source The Array2D instance to be copied.
       */
      Kernel(const brick::numeric::Array2D<ElementType>& source);


      /**
       * Construct a kernel around external data.  Kernels constructed in
       * this way will copy the data from the input pointer.
       *
       * @param rows Number of rows in the kernel after successful
       * construction.
       *
       * @param columns Number of columns in the kernel after successful
       * construction.
       *
       * @param dataPtr A C-style array of ElementType from which the newly
       * constructed Kernel should copy its data.
       */
      Kernel(size_t rows, size_t columns, ElementType* const dataPtr);


      /**
       * This constructor allows us to implicitly make a separable
       * Kernel instance from a pair of Array1D instances.  As with the
       * copy constructor, the newly created kernel has its own copy of
       * the data from the copied arrays.
       *
       * @param rowArray This argument specifies the portion of the
       * kernel which will be applied along image rows.
       *
       * @param columnArray This argument specifies the portion of the
       * kernel which will be applied along image columns.
       */
      Kernel(const brick::numeric::Array1D<ElementType>& rowArray,
             const brick::numeric::Array1D<ElementType>& columnArray);


      /**
       * Destroys the Kernel instance and deletes the internal data
       * store.
       */
      virtual
      ~Kernel();


      /**
       * The assignment operator does a deep copy.
       *
       * @param source The Kernel instance to be copied.
       */
      Kernel<ELEMENT_TYPE>&
      operator=(Kernel<ELEMENT_TYPE> const& source);


      /**
       * Return a copy of the kernel data in an Array2D object.
       *
       * @return The return value is an Array2D instance.  It does not
       * reference the same physical memory as the kernel.
       */
      brick::numeric::Array2D<ElementType>
      getArray2D() const;


      /**
       * This method is only valid for separable kernels; it returns
       * the separable kernel component which is parallel to the columns
       * of the image.
       *
       * @return The return value is an Array1D instance representing an
       * N-row by 1-column kernel component.
       */
      brick::numeric::Array1D<ElementType>
      getColumnComponent() const;


      /**
       * This member function returns the number of columns in the kernel.
       *
       * @return The return value is the width, in columns, of the kernel.
       */
      size_t getColumns() const;


      /**
       * This method is only valid for separable kernels; it returns
       * the separable kernel component which is parallel to the rows
       * of the image.
       *
       * @return The return value is an Array1D instance representing a
       * 1-row by M-column kernel component.
       */
      brick::numeric::Array1D<ElementType>
      getRowComponent() const;


      /**
       * This member function returns the number of rows in the kernel.
       *
       * @return The return value is the height, in rows, of the kernel.
       */
      size_t getRows() const;


      /**
       * This method returns true if the kernel is separable, false
       * otherwise.
       *
       * @return The return value indicates whether or no the kernel is
       * separable.
       */
      bool
      isSeparable() const {return m_isSeparable;}


      /**
       * This method allows the user to set the contents of the kernel.
       * This does a deep copy of the input array.
       *
       * @param dataArray This argument specifies the kernel data to be
       * copied.
       */
      void
      setArray2D(const brick::numeric::Array2D<ElementType>& dataArray);


      /**
       * This method allows the user to set the contents of a separable
       * kernel.  Calling this method forces the kernel instance to
       * become separable.  This method does a deep copy of the input
       * arrays.
       *
       * @param rowArray This argument specifies the separable kernel
       * component which is parallel to the rows of the image.  It
       * represents a 1-row by M-column kernel component.
       *
       * @param columnArray This argument specifies the separable kernel
       * component which is parallel to the columns of the image.  It
       * represents an N-row by 1-column kernel component.
       */
      void
      setSeparableComponents(
        const brick::numeric::Array1D<ElementType>& rowArray,
        const brick::numeric::Array1D<ElementType>& columnArray);


    private:

      brick::numeric::Array2D<ElementType> m_data;
      bool m_isSeparable;
      size_t m_separableColumns;
      size_t m_separableRows;

    };

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/kernel_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KERNEL_HH */
