/**
***************************************************************************
* @file brick/numeric/bSpline_impl.hh
*
* Header file defining inline and template functions from bSpline.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_BSPLINE_IMPL_HH
#define BRICK_NUMERIC_BSPLINE_IMPL_HH

// This file is included by bSpline.hh, and should not be directly included
// by user code, so no need to include bSpline.hh here.
//
// #include <brick/numeric/bSpline.hh>

#include <cmath>
#include <brick/numeric/functional.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace numeric {

    // Note(xxx): Remove m_numberOfNodes?

    // This constructor builds a BSpline instance of unspecified length.
    template <class Type>
    BSpline<Type>::
    BSpline(size_t order, bool isPeriodic)
      : m_coefficientMatrixArray(),
        m_controlPoints(),
        m_controlPointVectorArray(),
        m_cumulativeKnotCounts(),
        m_inputVector(order + 1),
        m_isPeriodic(isPeriodic),
        m_isUniform(true),
        m_knotPositionArray(),
        m_nodePositionArray(),
        m_numberOfNodes(0),
        m_order(order),
        m_orderPlusOne(order + 1)
    {
      // Empty
    }


    // The copy constructor does a deep copy.
    template <class Type>
    BSpline<Type>::
    BSpline(const BSpline& other)
      : m_coefficientMatrixArray(other.m_coefficientMatrixArray.copy()),
        m_controlPoints(other.m_controlPoints),
        m_controlPointVectorArray(other.m_controlPointVectorArray.copy()),
        m_cumulativeKnotCounts(other.m_cumulativeKnotCounts.copy()),
        m_inputVector(other.m_inputVector.copy()),
        m_isPeriodic(other.isPeriodic),
        m_isUniform(other.isUniform),
        m_knotPositionArray(other.m_knotPositionArray.copy()),
        m_nodePositionArray(other.m_nodePositionArray.copy()),
        m_numberOfNodes(other.m_numberOfNodes),
        m_order(other.m_order),
        m_orderPlusOne(other.m_orderPlusOne)
    {
      // Empty
    }




    // This member function returns the maximum value for the spline
    // parameter S.
    template <class Type>
    double
    BSpline<Type>::
    getMaximumSValue()
    {
      if(m_nodePositionArray.size() == 0) {
        BRICK_THROW(brick::common::StateException,
                    "BSpline::getMaximumSValue()",
                    "Node positions have not been set.  You must first call"
                    " setNumberOfNodes() or setNodePositions().");
      }
      return m_nodePositionArray[m_nodePositionArray.size() - 1];
    }


    // This member function returns the minimum value for the spline
    // parameter S.
    template <class Type>
    double
    BSpline<Type>::
    getMinimumSValue()
    {
      if(m_nodePositionArray.size() == 0) {
        BRICK_THROW(brick::common::StateException,
                    "BSpline::getMinimumSValue()",
                    "Node positions have not been set.  You must first call"
                    " setNumberOfNodes() or setNodePositions().");
      }
      return m_nodePositionArray[0];
    }


    // This member function sets the values of the control points of
    // the spline.
    template <class Type>
    void
    BSpline<Type>::
    setControlPoints(const std::vector<Type>& controlPoints)
    {
      if(m_isPeriodic) {
        if(controlPoints.size() != m_numberOfNodes - 1) {
          std::ostringstream message;
          message << "Argument controlPoints has " << controlPoints.size()
                  << " elements, but should have " << m_numberOfNodes - 1
                  << " elements for a periodic " << m_numberOfNodes
                  << " node spline because the first and last nodes overlap.";
          BRICK_THROW(brick::common::ValueException,
                      "BSpline::setControlPoints()",
                      message.str().c_str());
        }
      } else {
        if(controlPoints.size() != m_numberOfNodes) {
          std::ostringstream message;
          message << "Argument controlPoints has " << controlPoints.size()
                  << " elements, but should have " << m_numberOfNodes
                  << " elements for a non-periodic " << m_numberOfNodes
                  << " node spline.";
          BRICK_THROW(brick::common::ValueException,
                      "BSpline::setControlPoints()",
                      message.str().c_str());
        }
      }
      if(m_knotPositionArray.size() < m_nodePositionArray.size()) {
        BRICK_THROW(brick::common::StateException,
                    "BSpline::setControlPoints()",
                    "Knot positions have not been initialized.  You must first"
                    " call setKnotMultiplicies(), or else call"
                    " setNumberOfNodes() with argument setKnotMultipliciesFlag"
                    " set to true.");
      }
      m_controlPoints = controlPoints;
      m_controlPointVectorArray = this->getControlPointVectors(
        m_order, m_numberOfNodes, m_cumulativeKnotCounts, controlPoints);
    }


    // This member function sets the knot multiplicity at each node
    // of the spline.
    template <class Type>
    void
    BSpline<Type>::
    setKnotMultiplicities(const std::vector<size_t>& knotMultiplicities)
    {
      if(knotMultiplicities.size() != m_numberOfNodes) {
        std::ostringstream message;
        message << "Argument knotMultiplicities has "
                << knotMultiplicities.size() << " elements, when a(n) "
                << m_numberOfNodes << " element vector was expected.";
        BRICK_THROW(brick::common::ValueException,
                    "BSpline::setKnotMultiplicities()",
                    message.str().c_str());
      }
      if(m_isPeriodic) {
        if(knotMultiplicities[0] != knotMultiplicities[m_numberOfNodes - 1]) {
          std::ostringstream message;
          message << "For a periodic spline, first and last node are actually "
                  << "the same and so must have the same knot multiplicity, "
                  << "but knotMultiplicities[0] = " << knotMultiplicities[0]
                  << " and knotMultiplicities[" << m_numberOfNodes - 1
                  << "] = " << knotMultiplicities[m_numberOfNodes - 1] << ".";
          BRICK_THROW(brick::common::ValueException,
                      "BSpline::setKnotMultiplicities()",
                      message.str().c_str());
        }
      } else {
        if(knotMultiplicities[0] != m_order
           || knotMultiplicities[m_numberOfNodes - 1] != m_order) {
          std::ostringstream message;
          message << "For a non-periodic spline of order " << m_order << ", "
                  << "first and last knot must have multiplicity " << m_order
                  << ", but knotMultiplicities[0] = " << knotMultiplicities[0]
                  << " and knotMultiplicities[" << m_numberOfNodes - 1
                  << "] = " << knotMultiplicities[m_numberOfNodes - 1] << ".";
          BRICK_THROW(brick::common::ValueException,
                      "BSpline::setKnotMultiplicities()",
                      message.str().c_str());
        }
      }
      if(m_cumulativeKnotCounts.size() != knotMultiplicities.size()) {
        m_cumulativeKnotCounts.reinit(knotMultiplicities.size());
      }
      std::partial_sum(knotMultiplicities.begin(), knotMultiplicities.end(),
                       m_cumulativeKnotCounts.begin());
      size_t numberOfKnots = m_cumulativeKnotCounts[m_numberOfNodes - 1];
      if(m_knotPositionArray.size() != numberOfKnots) {
        m_knotPositionArray.reinit(numberOfKnots);
      }
      size_t knotIndex = 0;
      for(size_t nodeIndex = 0; nodeIndex < m_numberOfNodes; ++nodeIndex) {
        while(knotIndex < m_cumulativeKnotCounts[nodeIndex]) {
          m_knotPositionArray[knotIndex] = m_nodePositionArray[nodeIndex];
          ++knotIndex;
        }
      }
      m_coefficientMatrixArray = this->getCoefficientMatrices(
        m_order, m_numberOfNodes, m_nodePositionArray, m_knotPositionArray,
        m_cumulativeKnotCounts);
    }



    // This member function both specifies the number of nodes in
    // the spline and sets the node positions so that the spline is
    // "uniform".
    template <class Type>
    void
    BSpline<Type>::
    setNumberOfNodes(size_t numberOfNodes, bool setKnotMultiplicitiesFlag)
    {
      std::vector<double> nodePositions(numberOfNodes);
      for(size_t index0 = 0; index0 < numberOfNodes; ++index0) {
        nodePositions[index0] = static_cast<double>(index0);
      }
      this->setNumberOfNodes(
        numberOfNodes, nodePositions, setKnotMultiplicitiesFlag);
      m_isUniform = true;
    }


    // This member function specifies the number of nodes in the
    // spline and allows the user to set the position of each node.
    template <class Type>
    void
    BSpline<Type>::
    setNumberOfNodes(size_t numberOfNodes,
                     const std::vector<double>& nodePositions,
                     bool setKnotMultiplicitiesFlag)
    {
      m_numberOfNodes = numberOfNodes;
      m_isUniform = false;
      this->setNodePositions(nodePositions);
      if(setKnotMultiplicitiesFlag) {
        std::vector<size_t> knotMultiplicities(numberOfNodes, size_t(1));
        if(!m_isPeriodic) {
          knotMultiplicities[0] = m_order + 1;
          knotMultiplicities[numberOfNodes - 1] = m_order + 1;
        }
        this->setKnotMultiplicities(knotMultiplicities);
      }
    }


    // The assigment operator does a deep copy.
    template <class Type>
    BSpline<Type>&
    BSpline<Type>::
    operator=(const BSpline<Type>& other)
    {
      if(&other != this) {
        m_coefficientMatrixArray = other.m_coefficientMatrixArray.copy();
        m_controlPoints = other.m_controlPoints;
        m_controlPointVectorArray = other.m_controlPointVectorArray.copy();
        m_cumulativeKnotCounts = other.m_cumulativeKnotCounts.copy();
        m_inputVector = other.m_inputVector.copy();
        m_isPeriodic = other.isPeriodic;
        m_isUniform = other.isUniform;
        m_knotPositionArray = other.m_knotPositionArray.copy();
        m_nodePositionArray = other.m_nodePositionArray.copy();
        m_numberOfNodes = other.m_numberOfNodes;
        m_order = other.m_order;
        m_orderPlusOne = other.m_orderPlusOn;
      }
    }


    // This operator evaluates the spline at the specified value of
    // spline parameter s.
    template <class Type>
    Type
    BSpline<Type>::
    operator()(double sValue)
    {
      // Construct a vector of powers of s.
      double accumulator = 1.0;
      m_inputVector[0] = accumulator;
      for(size_t index0 = 1; index0 < m_inputVector.size(); ++index0) {
        accumulator *= sValue;
        m_inputVector[index0] = accumulator;
      }

      // Select the appropriate coefficients and control points for
      // the span to which sValue belongs.
      size_t spanNumber = this->getSpanNumber(sValue);
      const Array1D<Type>& controlPointVector =
        m_controlPointVectorArray[spanNumber];
      const Array2D<double>& coefficientMatrix =
        m_coefficientMatrixArray[spanNumber];


      // Actually compute the polynomial values and multiply by
      // control points.
      Type result = (dot<Type>(coefficientMatrix.row(0), m_inputVector)
                     * controlPointVector[0]);
      for(size_t index1 = 1; index1 < m_orderPlusOne; ++index1) {
        result += (dot<Type>(coefficientMatrix.row(index1), m_inputVector)
                   * controlPointVector[index1]);
      }
      return result;
    }


    // This protected member function computes one spline basis
    // function for use in calculating spline values.
    template <class Type>
    Polynomial<double>
    BSpline<Type>::
    computeBasisFunction(size_t order,
                         size_t spanNumber,
                         size_t componentNumber,
                         const Array1D<size_t>& cumulativeKnotCounts,
                         const Array1D<double>& knotPositions)
    {
      // Initialize 0-order basis functions according to the rule:
      //
      //   B_(n,1)(s) = 1     if k_n <= s < k_(n+1)
      //   B_(n,1)(s) = 0     otherwise
      //
      // Where B_(n,1)(s) is the 0-order basis function with support
      // (k_n <= s < k_(n+1), k_n is the position of the n^th knot in
      // the spline, k_(n+1) is the position of the (n+1)^st knot in
      // the spline, and s is the parametric distance along the
      // spline.
      //
      // Note that only one of the 0 order components is non-zero at
      // any given point along the spline.  The next three lines
      // create an array of order + 1 polynomials, initialize all of
      // them to zero, and then set the appropriate polynomial to one.
      Array1D< Polynomial<double> > componentArray(order + 1);
      componentArray = Polynomial<double>(0.0);
      componentArray[order - componentNumber] = Polynomial<double>(1.0);

      // Recursively compute higher order polynomials accoring to the rule:
      //
      //  B_(n,d)(s) = ((s - k_n)B_(n,d-1)(s)/(k_(n+d-1) - k_n)
      //                + (k_(n+d) - s)B_(n+1,d-1)(s)/(k_(n+d) - k_(n+1)))
      //
      // For each stage of the recursion.
      for(int orderIndex = 1; orderIndex <= static_cast<int>(order); ++orderIndex) {

        // Polynomials of the next higher order will be put into
        // newComponentArray.
        Array1D< Polynomial<double> > newComponentArray(
          componentArray.size() - 1);

        // For each element of newComponentArray.
        for(size_t subComponentIndex = 0;
            subComponentIndex < newComponentArray.size();
            ++subComponentIndex) {
          // knotNumber is the number of the first (lowest numbered)
          // knot in the recursion rule.  knotNumber corresponds to
          // "n" in the equation above.
          int knotNumber =
            (static_cast<int>(cumulativeKnotCounts[spanNumber]) - 1
             - static_cast<int>(m_order)
             + static_cast<int>(componentNumber)
             + static_cast<int>(subComponentIndex));
          double k_n =
            this->getKnotPosition(knotNumber, knotPositions);
          double k_nPlus1 =
            this->getKnotPosition(knotNumber + 1, knotPositions);
          double k_nPlusDMinus1 =
            this->getKnotPosition(knotNumber + orderIndex, knotPositions);
          double k_nPlusD =
            this->getKnotPosition(knotNumber + orderIndex + 1, knotPositions);
          double span0 = k_nPlusDMinus1 - k_n;
          double span1 = k_nPlusD - k_nPlus1;
          Polynomial<double> scaleFunction0(1.0 / span0, -k_n / span0);
          Polynomial<double> scaleFunction1(-1.0 / span1, k_nPlusD / span1);
          newComponentArray[subComponentIndex] = (
            scaleFunction0 * componentArray[subComponentIndex]
            + scaleFunction1 * componentArray[subComponentIndex + 1]);
        }
        // Now that the recursive rule has been carried out for this
        // level, (this value of "d" in the equation above), we can
        // forget the previous level and get ready to recurse.
        componentArray = newComponentArray;
      }
      if(componentArray.size() != 1) {
        BRICK_THROW(brick::common::LogicException,
                    "BSpline::computeBasisFunction()",
                    "Recursion terminated incorrectly.");
      }
      return componentArray[0];
    }


    // This protected member function returns an array in which each
    // element corresponds to one span of the spline, and contains a
    // matrix of the polynomial coefficients of the basis functions
    // that affect the spline values within that span.
    template <class Type>
    Array1D< Array2D<double> >
    BSpline<Type>::
    getCoefficientMatrices(size_t order,
                           size_t numberOfNodes,
                           const Array1D<double>& /* nodePositions */,
                           const Array1D<double>& knotPositions,
                           const Array1D<size_t>& cumulativeKnotCounts)
    {
      size_t numberOfSpans = numberOfNodes - 1;
      Array1D< Array2D<double> > coefficientMatrixArray(numberOfSpans);
      for(size_t spanNumber = 0; spanNumber < numberOfSpans; ++spanNumber) {
        coefficientMatrixArray[spanNumber] =
          Array2D<double>(order + 1, order + 1);

        for(size_t componentIndex = 0; componentIndex < order + 1;
            ++componentIndex) {
          Polynomial<double> basisFunction = this->computeBasisFunction(
            order, spanNumber, componentIndex, cumulativeKnotCounts,
            knotPositions);
          if(basisFunction.getOrder() != order) {
            BRICK_THROW(brick::common::LogicException,
                        "BSpline::getCoefficientMatrices()",
                        "Basis function has incorrect order.");
          }
          coefficientMatrixArray[spanNumber].row(componentIndex).copy(
            basisFunction.getCoefficientArray());
        }
      }
      return coefficientMatrixArray;
    }


    // This protected member function returns an array in which each
    // element corresponds to one span of the spline, and contains
    // the control-point values that affect the spline values
    // within that span.  This function is used to allow efficient
    template <class Type>
    Array1D< Array1D<double> >
    BSpline<Type>::
    getControlPointVectors(size_t order, size_t numberOfNodes,
                           const Array1D<size_t>& cumulativeKnotCounts,
                           const std::vector<Type>& controlPoints)
    {
      size_t numberOfSpans = numberOfNodes - 1;
      size_t numberOfBasisFunctions = order + 1;
      size_t numberOfKnots;
      if(m_isPeriodic) {
        numberOfKnots = cumulativeKnotCounts[numberOfNodes - 2];
      } else {
        numberOfKnots = cumulativeKnotCounts[numberOfNodes - 1];
      }
      Array1D< Array1D<Type> > controlPointVectorArray(numberOfSpans);
      for(size_t spanNumber = 0; spanNumber < numberOfSpans; ++spanNumber) {
        Array1D<Type> newVector(numberOfBasisFunctions);
        int knotIndex =
          static_cast<int>(cumulativeKnotCounts[spanNumber])
		  - static_cast<int>(numberOfBasisFunctions);
        for(size_t controlPointIndex = 0;
            controlPointIndex < newVector.size();
            ++controlPointIndex) {
          if(knotIndex < 0) {
            newVector[controlPointIndex] =
              controlPoints[numberOfKnots + knotIndex];
          } else if(knotIndex >= static_cast<int>(numberOfKnots)) {
            newVector[controlPointIndex] =
              controlPoints[knotIndex - numberOfKnots];
          } else {
            newVector[controlPointIndex] = controlPoints[knotIndex];
          }
          ++knotIndex;
        }
        controlPointVectorArray[spanNumber] = newVector;
      }
      return controlPointVectorArray;
    }


    // This protected member function wraps argument knotNumber so
    // that it is in the range [0, knotPositions.size() - 1], and then
    // returns the corresponding value from knotPositions.
    template <class Type>
    double
    BSpline<Type>::
    getKnotPosition(int knotNumber,
                    const Array1D<double>& knotPositions)
    {
      if(knotNumber < 0) {
        int wrappedIndex =
          static_cast<int>(knotPositions.size()) - 1 + knotNumber;
        double offsetDistance = (knotPositions[knotPositions.size() - 1]
                                 - knotPositions[wrappedIndex]);
        return knotPositions[0] - offsetDistance;
      } else if(knotNumber >= static_cast<int>(knotPositions.size())) {
        int wrappedIndex =
          knotNumber - static_cast<int>(knotPositions.size()) + 1;
        double offsetDistance = (knotPositions[wrappedIndex]
                                 - knotPositions[0]);
        return knotPositions[knotPositions.size() - 1] + offsetDistance;
      }
      return knotPositions[knotNumber];
    }


    // This protected member function returns the number of the span
    // in which the specified spline parameter value lies.
    template <class Type>
    size_t
    BSpline<Type>::
    getSpanNumber(double sValue)
    {
      size_t returnValue;
      if(m_isUniform) {
        int maximumSpanNumber = static_cast<int>(m_numberOfNodes) - 1;
        int spanNumber = int(floor(sValue));
        if(m_isPeriodic) {
          while(spanNumber < 0) {
            spanNumber += maximumSpanNumber;
          }
          while(spanNumber >= maximumSpanNumber) {
            spanNumber -= maximumSpanNumber;
          }
        } else {
          if((spanNumber < 0) || (spanNumber >= maximumSpanNumber)) {
            BRICK_THROW(brick::common::ValueException,
                        "BSpline::getSpanNumber()",
                        "Argument sValue is out of range.");
          }
        }
        returnValue = static_cast<size_t>(spanNumber);
      } else { // if(m_isUniform)
        double maximumSValue = m_nodePositionArray[m_numberOfNodes - 1];
        double minimumSValue = m_nodePositionArray[0];
        if(m_isPeriodic) {
          while(sValue < minimumSValue) {
            sValue += (maximumSValue - minimumSValue);
          }
          while(sValue >= maximumSValue) {
            sValue -= (maximumSValue - minimumSValue);
          }
        } else {
          if((sValue < minimumSValue) || (sValue >= maximumSValue)) {
            BRICK_THROW(brick::common::ValueException,
                        "BSpline::getSpanNumber()",
                        "Argument sValue is out of range.");
          }
        }
        Array1D<double>::const_iterator nodeIterator = std::upper_bound(
          m_nodePositionArray.begin(), m_nodePositionArray.end(), sValue);
        returnValue =
          static_cast<size_t>(nodeIterator - m_nodePositionArray.begin());
      }
      return returnValue;
    }


    // This protected member function sets the positions of the
    // nodes in the spline.
    template <class Type>
    void
    BSpline<Type>::
    setNodePositions(const std::vector<double>& nodePositions)
    {
      if(nodePositions.size() != m_numberOfNodes) {
        std::ostringstream message;
        message << "Argument nodePositions has "
                << nodePositions.size() << " elements, when a(n) "
                << m_numberOfNodes << " element vector was expected.";
        BRICK_THROW(brick::common::ValueException,
                    "BSpline::setNodePositions()",
                    message.str().c_str());
      }
      if(m_nodePositionArray.size() != nodePositions.size()) {
        m_nodePositionArray.reinit(nodePositions.size());
      }
      std::copy(nodePositions.begin(), nodePositions.end(),
                m_nodePositionArray.begin());
    }

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_BSPLINE_IMPL_HH */
