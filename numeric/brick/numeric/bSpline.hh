/**
***************************************************************************
* @file brick/numeric/bSpline.hh
*
* Header file declaring the BSpline class.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_BSPLINE_HH
#define BRICK_NUMERIC_BSPLINE_HH

#include <vector>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/polynomial.hh>

namespace brick {

  namespace numeric {

    /**
     ** Warning: The test suite for this class is still incomplete.
     ** It may contain bugs.  In particular, the non-periodic
     ** implementation is known to be broken.
     **
     ** This class template implements a B-spline of arbitrary order.
     ** It is templated on control point type so that you can, for
     ** example, create a spline with 2D or 3D values by by specifying
     ** control points of type Vector2D or Vector3D.  Note that the
     ** splines represented by this class are all 1D in the sense that
     ** you can represent them using a parametric function of one
     ** parameter.  If you set the template parameter to be Vector2D,
     ** then you wind up with a function that converts this 1D input
     ** parameter into 2D output values.  If you want a spline that
     ** has two input parameters, perhaps to represent a surface, then
     ** you need to use the similar, but more restrictive, BSpline2D
     ** class.  The BSpline class supports both periodic and
     ** non-periodic splines, supports both uniform and non-uniform
     ** node spacing, and supports multiple knots at each node to
     ** allow corners and discontinuities.
     **/
    template <class Type>
    class BSpline {
    public:

      /** 
       * This constructor builds a BSpline instance of unspecified length.
       * 
       * @param order This argument sets the order of the spline.  For
       * a quadratic spline, set order to 2.  For a cubic spline, set
       * order to 3.
       * 
       * @param isPeriodic If this argument is true, the spline will
       * be periodic.  That is, its last node will overlap its first,
       * and the spline parameter will wrap around from its maximum
       * value back to zero.
       */
      BSpline(size_t order=2, bool isPeriodic=true);


      /** 
       * The copy constructor does a deep copy.
       * 
       * @param other This argument is the BSpline instance to be copied.
       */
      BSpline(const BSpline& other);


      /** 
       * This member function returns the maximum value for the spline
       * parameter S.  For a non-periodic spline, calling
       * operator()(double) with an argument greater than or equal to
       * the result of getMaximumSValue() is an error.  For a periodic
       * spline, calling operator()(double) with an argument greater
       * than or equal to the result of getMaximumSValue() is not an
       * error, but the parameter will be wrapped around to the
       * beginning of the spline.
       * 
       * @return The return value is the position of the last node in
       * the spline.
       */
      double
      getMaximumSValue();


      /** 
       * This member function returns the minimum value for the spline
       * parameter S.  For a non-periodic spline, calling
       * operator()(double) with an argument less than the result of
       * getMinimumSValue() is an error.  For a periodic spline,
       * calling operator()(double) with an argument less than the
       * result of getMinimumSValue() is not an error, but the
       * parameter will be wrapped around to the end of the spline.
       * 
       * @return The return value is the position of the first node in
       * the spline.
       */
      double
      getMinimumSValue();
      

      /** 
       * This member function sets the values of the control points of
       * the spline.  If the spline is periodic, then the value of the
       * final control point should be omitted; it will be
       * automatically copied from the value of the first control
       * point.
       * 
       * @param controlPoints This argument specifies the control
       * point values for the spline.
       */
      void
      setControlPoints(const std::vector<Type>& controlPoints);


      /** 
       * This member function sets the knot multiplicity at each node
       * of the spline.  Setting the knot multiplicity of a node to N
       * will introduce, at that node, a discontinuity in the ((order
       * - N) + 1)th derivative of the spline.  This means that a knot
       * multipicity equal to the order of the spline will introduce a
       * "corner" (discontinuity in the 1st derivative) and a
       * multipicity equal to (order + 1) will introduce a break in
       * the spline.  By default, all knot multiplicities are set to
       * 1, except for the first and last nodes of a non-periodic
       * spline, which are set to (order + 1).
       * 
       * @param knotMultiplicities This argument specifies the knot
       * multipicity at each node in the spline.  For a periodic
       * spline, the first and last nodes overlap, and must have the
       * same knot multiplicity.  For a non-periodic spline, the first and
       * last nodes must have knot multiplicity equal to (order + 1).
       */
      void
      setKnotMultiplicities(const std::vector<size_t>& knotMultiplicities);


      /** 
       * This member function both specifies the number of nodes in
       * the spline and sets the node positions so that the spline is
       * "uniform".  The node positions will be set so that the first
       * node lies at spline parameter s = 0.0, the second node lies
       * at s = 1.0, the third at s = 2.0, and so on.  For periodic
       * splines, the first and last node overlap, and represent the
       * same control point.  For both periodic and non-periodic
       * splines, the number of spans will be equal to numberOfNodes -
       * 1.
       * 
       * @param numberOfNodes This argument specifies how many nodes
       * the spline should have.  For periodic splines, the first and
       * last nodes represent the same physical point on the spline.
       * 
       * @param setKnotMultiplicitiesFlag This argument is used to
       * avoid redundant calculations if the knot multiplicities will
       * be explicitly set later.  If this argument is set to true,
       * the knot multiplicity of each node will be set to the default
       * value.  The default value is 1 for each node except the first
       * and last nodes of a non-periodic spline, for which the
       * default knot multiplicity is (order + 1).  If this argument
       * is set to false, the knot multiplicities will not be set, and
       * the calling context must explicitly set them by calling
       * member function setKnotMultiplicities().
       */
      void
      setNumberOfNodes(size_t numberOfNodes,
                       bool setKnotMultiplicitiesFlag=true);
      

      /** 
       * This member function specifies the number of nodes in the
       * spline and allows the user to set the position of each node.
       * The node position values must be monotonically increasing
       * with node number.  For periodic splines, the first and last
       * node overlap, and represent the same control point, however
       * different values should be specified for the positions of the
       * first and last nodes.  When the spline parameter, s, reaches
       * the position of the final node, it will be wrapped around, as
       * if it were actually set to the position of the first node.
       * For both periodic and non-periodic splines, the number of
       * spans will be equal to numberOfNodes - 1.
       * 
       * @param numberOfNodes This argument specifies how many nodes
       * the spline should have.  For periodic splines, the first and
       * last nodes represent the same physical point on the spline.
       * 
       * @param nodePositions This argument specifies the positions
       * (in spline parameter space) of the nodes.  This argument can
       * be used to create non-uniform splines.
       * 
       * @param setKnotMultiplicitiesFlag This argument is used to
       * avoid redundant calculations if the knot multiplicities will
       * be explicitly set later.  If this argument is set to true,
       * the knot multiplicity of each node will be set to the default
       * value.  The default value is 1 for each node except the first
       * and last nodes of a non-periodic spline, for which the
       * default knot multiplicity is (order + 1).  If this argument
       * is set to false, the knot multiplicities will not be set, and
       * the calling context must explicitly set them by calling
       * member function setKnotMultiplicities().
       */
      void
      setNumberOfNodes(size_t numberOfNodes,
                       const std::vector<double>& nodePositions,
                       bool setKnotMultiplicitiesFlag=true);
      

      
      /** 
       * The assigment operator does a deep copy.
       * 
       * @param other This argument is the BSpline instance to be copied.
       */
      BSpline<Type>&
      operator=(const BSpline<Type>& other);

      
      /** 
       * This operator evaluates the spline at the specified value of
       * spline parameter s.
       * 
       * @return The return value is the calculated spline value.
       */
      Type
      operator()(double sValue);
      

    protected:

      /** 
       * This protected member function computes one spline basis
       * function for use in calculating spline values.
       * 
       * @param order This argument specifies the order of the spline
       * (2 for quadratic, 3 for cubic, etc.)
       * 
       * @param spanNumber This argument specifies the span for which
       * the basis function is being computed.
       * 
       * @param componentNumber This argument specifies which of the
       * (order + 1) basis functions that overlap the span is to be
       * computed.
       * 
       * @param cumulativeKnotCounts This argument specifies, for each
       * node, the total number of knots at that node plus the total
       * number of knots at preceding nodes.
       * 
       * @param knotPositions This argument specifies the position
       * (spline parameter) of each knot.  If the spline contains
       * nodes with multiple knots, then the argument will contain
       * consecutive entries with the same value.
       * 
       * @return The return value is the requested polynomial.
       */
      Polynomial<double>
      computeBasisFunction(size_t order,
                           size_t spanNumber,
                           size_t componentNumber,
                           const Array1D<size_t>& cumulativeKnotCounts,
                           const Array1D<double>& knotPositions);

      
      /** 
       * This protected member function returns an array in which each
       * element corresponds to one span of the spline, and contains
       * the control-point values that affect the spline values
       * within that span.  This function is used to allow efficient
       * calculation of spline values.
       * 
       * @param order This argument is the order of the spline.
       * 
       * @param numberOfNodes This argument specifies the number of
       * nodes in the spline.
       * 
       * @param cumulativeKnotCounts This argument specifies, for each
       * node, the total number of knots at that node plus the total
       * number of knots at preceding nodes.
       * 
       * @param controlPoints This argument specifies the actual
       * control point values.
       * 
       * @return The return value is an array of arrays of
       * pre-selected control point values.
       */
      Array1D< Array1D<double> >
      getControlPointVectors(size_t order,
                             size_t numberOfNodes,
                             const Array1D<size_t>& cumulativeKnotCounts,
                             const std::vector<Type>& controlPoints);

      
      /** 
       * This protected member function returns an array in which each
       * element corresponds to one span of the spline, and contains a
       * matrix of the polynomial coefficients of the basis functions
       * that affect the spline values within that span.  This
       * function is used to allow efficient calculation of spline
       * values.
       * 
       * @param order This argument is the order of the spline.
       * 
       * @param numberOfNodes This argument specifies the number of
       * nodes in the spline.
       *
       * @param nodePositions This argument specifies the positions of
       * each node in the spline.  Positions are defined in terms of
       * spline parameter s.
       *
       * @param knotPositions This argument specifies the positions of
       * each knot in the spline.  Each knot has the same position as
       * the node with which it is associated.  Positions are defined
       * in terms of spline parameter s.
       *
       * @param cumulativeKnotCounts This argument specifies, for each
       * node, the total number of knots at that node plus the total
       * number of knots at preceding nodes.
       * 
       * @return The return value is an array of 2D arrays of
       * coefficient values.
       */
      Array1D< Array2D<double> >
      getCoefficientMatrices(size_t order,
                             size_t numberOfNodes,
                             const Array1D<double>& nodePositions,
                             const Array1D<double>& knotPositions,
                             const Array1D<size_t>& cumulativeKnotCounts);


      /** 
       * This protected member function wraps argument knotNumber so
       * that it is in the range [0, knotPositions.size() - 1], and
       * then returns the corresponding value from knotPositions.
       * 
       * @param knotNumber This argument is the number of the knot,
       * possibly out-of-range and needing to be wrapped.
       * 
       * @param knotPositions This argument is an array of knot positions.
       * 
       * @return The return value is the appropriate value from
       * argument knotPositions.
       */
      double
      getKnotPosition(int knotNumber,
                      const Array1D<double>& knotPositions);


      /** 
       * This protected member function returns the number of the span
       * in which the specified spline parameter value lies.
       * 
       * @param sValue This argument indicates the point of interest
       * along the spline.
       * 
       * @return The return value is the corresponding span number.
       */
      size_t
      getSpanNumber(double sValue);
      

      /** 
       * This protected member function sets the positions of the
       * nodes in the spline.
       * 
       * @param nodePositions This argument is a vector of node
       * positions.
       */
      void
      setNodePositions(const std::vector<double>& nodePositions);


      Array1D< Array2D<double> > m_coefficientMatrixArray;
      std::vector<Type> m_controlPoints;
      Array1D< Array1D<Type> > m_controlPointVectorArray;
      Array1D<size_t> m_cumulativeKnotCounts;
      Array1D<double> m_inputVector;
      bool m_isPeriodic;
      bool m_isUniform;
      Array1D<double> m_knotPositionArray;
      Array1D<double> m_nodePositionArray;
      size_t m_numberOfNodes;
      size_t m_order;
      size_t m_orderPlusOne;
    };


  } // namespace numeric
  
} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/bSpline_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_BSPLINE_HH */
