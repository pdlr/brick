/**
***************************************************************************
* @file brick/computerVision/kernel_impl.hh
*
* Header file defining Kernel class inline and template members.
*
* Copyright (C) 2006, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KERNEL_IMPL_HH
#define BRICK_COMPUTERVISION_KERNEL_IMPL_HH

// This file is included by kernel.hh, and should not be directly included
// by user code, so no need to include kernel.hh here.
//
// #include <brick/computerVision/kernel.hh>


namespace brick {

  namespace computerVision {

    // Default constructor initializes to zero size.
    template <class ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>::
    Kernel()
      : m_data(),
        m_isSeparable(false),
        m_separableColumns(0),
        m_separableRows(0)
    {
      // Empty.
    }


    // The copy constructor does a deep copy.
    template <class ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>::
    Kernel(const Kernel<ELEMENT_TYPE> &source)
      : m_data(source.m_data.copy()),
        m_isSeparable(source.m_isSeparable),
        m_separableColumns(source.m_separableColumns),
        m_separableRows(source.m_separableRows)
    {
      // Empty.
    }


    // This constructor allows us to implicitly make a Kernel instance
    // from an Array2D using a deep copy.
    template <class ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>::
    Kernel(const brick::numeric::Array2D<ElementType> &source)
      : m_data(source.copy()),
        m_isSeparable(false),
        m_separableColumns(0),
        m_separableRows(0)
    {
      // Empty.
    }


    // Construct a kernel around external data using a deep copy.
    template <class ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>::
    Kernel(size_t rows, size_t columns, ElementType* const dataPtr)
      : m_data(rows, columns),
        m_isSeparable(false),
        m_separableColumns(0),
        m_separableRows(0)
    {
      m_data.copy(dataPtr);
    }


    // This constructor allows us to implicitly make a separable
    // Kernel instance from a pair of Array1D instances.
    template <class ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>::
    Kernel(const brick::numeric::Array1D<ElementType>& rowArray,
           const brick::numeric::Array1D<ElementType>& columnArray)
      : m_data(1, rowArray.size() + columnArray.size()),
        m_isSeparable(true),
        m_separableColumns(rowArray.size()),
        m_separableRows(columnArray.size())
    {
      std::copy(rowArray.begin(), rowArray.end(), m_data.begin());
      std::copy(columnArray.begin(), columnArray.end(),
                m_data.begin() + rowArray.size());
    }


    // Destroys the Kernel instance and deletes the internal data
    // store.
    template <class ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>::
    ~Kernel()
    {
      // Empty.
    }


    // The copy assignment operator does a deep copy.
    template <class ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>&
    Kernel<ELEMENT_TYPE>::
    operator=(const Kernel<ELEMENT_TYPE> &source)
    {
        m_data = source.m_data.copy();
        m_isSeparable = source.m_isSeparable;
        m_separableColumns = source.m_separableColumns;
        m_separableRows = source.m_separableRows;
        return *this;
    }


    // Return a copy of the kernel data in an Array2D object.
    template <class ELEMENT_TYPE>
    brick::numeric::Array2D<ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>::
    getArray2D() const
    {
      if(m_isSeparable) {
        brick::numeric::Array2D<ELEMENT_TYPE> synthesizedArray(
          m_separableRows, m_separableColumns);
        typename brick::numeric::Array2D<ELEMENT_TYPE>::iterator
          resultIterator = synthesizedArray.begin();
        for(size_t rowIndex = 0; rowIndex < m_separableRows; ++rowIndex) {
          for(size_t columnIndex = 0; columnIndex < m_separableColumns;
              ++columnIndex) {
            *resultIterator =
              m_data[columnIndex] * m_data[m_separableColumns + rowIndex];
	    ++resultIterator;
          }
        }
        return synthesizedArray;
      }
      // Not separable.
      return m_data.copy();
    }


    // This method is only valid for separable kernels; it returns
    // the separable kernel component which is parallel to the columns
    // of the image.
    template <class ELEMENT_TYPE>
    brick::numeric::Array1D<ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>::
    getColumnComponent() const
    {
      brick::numeric::Array1D<ELEMENT_TYPE> columnComponent(m_separableRows);
      columnComponent.copy(m_data.begin() + m_separableColumns);
      return columnComponent;
    }


    // This member function returns the number of columns in the kernel.
    template <class ELEMENT_TYPE>
    size_t
    Kernel<ELEMENT_TYPE>::
    getColumns() const
    {
      if(this->isSeparable()) {
        return m_separableColumns;
      }
      return m_data.columns();
    }


    // This method is only valid for separable kernels; it returns
    // the separable kernel component which is parallel to the rows
    // of the image.
    template <class ELEMENT_TYPE>
    brick::numeric::Array1D<ELEMENT_TYPE>
    Kernel<ELEMENT_TYPE>::
    getRowComponent() const
    {
      brick::numeric::Array1D<ELEMENT_TYPE> rowComponent(m_separableColumns);
      rowComponent.copy(m_data.begin());
      return rowComponent;
    }


    // This member function returns the number of rows in the kernel.
    template <class ELEMENT_TYPE>
    size_t
    Kernel<ELEMENT_TYPE>::
    getRows() const
    {
      if(this->isSeparable()) {
        return m_separableRows;
      }
      return m_data.rows();
    }


    // This method allows the user to set the contents of the kernel.
    template <class ELEMENT_TYPE>
    void
    Kernel<ELEMENT_TYPE>::
    setArray2D(const brick::numeric::Array2D<ElementType>& dataArray)
    {
      m_isSeparable = false;
      m_data = dataArray.copy();
      m_separableColumns = 0;
      m_separableRows = 0;
    }


    // This method allows the user to set the contents of a separable
    // kernel.
    template <class ELEMENT_TYPE>
    void
    Kernel<ELEMENT_TYPE>::
    setSeparableComponents(
      const brick::numeric::Array1D<ElementType>& rowArray,
      const brick::numeric::Array1D<ElementType>& columnArray)
    {
      m_isSeparable = true;
      m_data.reinit(1, rowArray.size() + columnArray.size());
      std::copy(rowArray.begin(), rowArray.end(), m_data.begin());
      std::copy(columnArray.begin(), columnArray.end(),
		m_data.begin() + rowArray.size());
      m_separableColumns = rowArray.size();
      m_separableRows = columnArray.size();
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KERNEL_IMPL_HH */
