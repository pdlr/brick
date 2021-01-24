/**
***************************************************************************
* @file brick/geometry/plane3D_impl.hh
*
* Source file defining the Plane3D class template.
*
* Copyright (C) 2007 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_PLANE3D_IMPL_HH
#define BRICK_GEOMETRY_PLANE3D_IMPL_HH

// This file is included by plane3D.hh, and should not be directly included
// by user code, so no need to include plane3D.hh here.
//
// #include <brick/geometry/plane3D.hh>

#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/subArray1D.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace geometry {

    // The default constructor initializes to the X-Y plane.
    template <class Type>
    Plane3D<Type>::
    Plane3D()
	: m_origin(0.0, 0.0, 0.0), m_directionVector0(1.0, 0.0, 0.0),
          m_directionVector1(0.0, 1.0, 0.0)
    {
      // Empty.
    }


    // This constructor initializes the plane using three points.
    template <class Type>
    Plane3D<Type>::
    Plane3D(const brick::numeric::Vector3D<Type>& point0,
            const brick::numeric::Vector3D<Type>& point1,
            const brick::numeric::Vector3D<Type>& point2,
            bool orthonormalize)
      : m_origin(point0),
        m_directionVector0(point1 - point0),
        m_directionVector1(point2 - point0)
    {
      if(orthonormalize) {
        m_directionVector0 /=
          brick::numeric::magnitude<Type>(m_directionVector0);
        m_directionVector1 -=
          (brick::numeric::dot<Type>(m_directionVector1, m_directionVector0)
           * m_directionVector0);
        m_directionVector1 /=
          brick::numeric::magnitude<Type>(m_directionVector1);
      }
    }


    // This constructor initializes the plane using a collection of
    // points.
    template <class Type>
    template <class Iterator>
    Plane3D<Type>::
    Plane3D(Iterator beginIterator,
            Iterator endIterator,
            double inlierPercentage)
      : m_origin(0.0, 0.0, 0.0),
        m_directionVector0(1.0, 0.0, 0.0),
        m_directionVector1(0.0, 1.0, 0.0)
    {
      // Note(xxx): Hasty implementation!  Clean up soon.
      std::vector< brick::numeric::Vector3D<Type> > allPointsVector;
      std::copy(beginIterator, endIterator,
                std::back_inserter(allPointsVector));

      int numberToInclude =
        static_cast<int>(allPointsVector.size() * inlierPercentage + 0.5);
      int numberToIgnore = allPointsVector.size() - numberToInclude;

      brick::numeric::Array1D<double> distances(allPointsVector.size());
      brick::numeric::Array1D<size_t> ignoredPoints(numberToIgnore);
      brick::numeric::Array1D< brick::numeric::Vector3D<Type> > currentPoints(
        allPointsVector.size());
      std::copy(allPointsVector.begin(), allPointsVector.end(),
                currentPoints.begin());
      bool isDone = false;
      ignoredPoints = 0;
      // xxx
      size_t iterationCount = 0;
      while(!isDone) {
        this->estimateFromSequence(currentPoints.begin(), currentPoints.end());
        for(size_t index0 = 0; index0 < allPointsVector.size(); ++index0) {
          distances[index0] = this->findDistance(allPointsVector[index0]);
        }
        brick::numeric::Array1D<size_t> pointIndices =
          brick::numeric::argsort(distances);
        if(static_cast<int>(currentPoints.size()) != numberToInclude) {
          currentPoints.reinit(numberToInclude);
        }
        for(int index1 = 0; index1 < numberToInclude; ++index1) {
          currentPoints[index1] = allPointsVector[pointIndices[index1]];
        }
        int index3 = 0;
        isDone = true;
        for(int index2 = numberToInclude;
            index2 < static_cast<int>(allPointsVector.size());
            ++index2) {
          if(ignoredPoints[index3] != pointIndices[index2]) {
            isDone = false;
            break;
          }
          ++index3;
        }
        // xxx
        if(++iterationCount >= 3) {
          break;
        }
        if(numberToInclude < static_cast<int>(allPointsVector.size())) {
          ignoredPoints = subArray(
            pointIndices, brick::numeric::Slice(numberToInclude, 0));
        } else {
          ignoredPoints.clear();
        }
      }
    }


    // The copy constructor deep copies its argument.
    template <class Type>
    Plane3D<Type>::
    Plane3D(const Plane3D<Type>& source)
      : m_origin(source.m_origin),
        m_directionVector0(source.m_directionVector0),
        m_directionVector1(source.m_directionVector1)
    {
      // Empty.
    }


    // Destructor.
    template <class Type>
    Plane3D<Type>::
    ~Plane3D()
    {
      // Empty.
    }


    // The assignment operator deep copies its argument.
    template <class Type>
    Plane3D<Type>&
    Plane3D<Type>::
    operator=(const Plane3D<Type>& source)
    {
      if(&source != this) {
        m_origin = source.m_origin;
        m_directionVector0 = source.m_directionVector0;
        m_directionVector1 = source.m_directionVector1;
      }
      return *this;
    }


    // This member function returns one of a pair of orthonormal
    // direction vectors that span the plane.
    template <class Type>
    brick::numeric::Vector3D<Type> const&
    Plane3D<Type>::
    getDirectionVector0() const
    {
      return m_directionVector0;
    }


    // This member function returns one of a pair of orthonormal
    // direction vectors that span the plane.
    template <class Type>
    brick::numeric::Vector3D<Type> const&
    Plane3D<Type>::
    getDirectionVector1() const
    {
      return m_directionVector1;
    }


    // xxx
    template <class Type>
    brick::numeric::Vector3D<Type>
    Plane3D<Type>::
    getNormal() const
    {
      return brick::numeric::cross(m_directionVector0, m_directionVector1);
    }


    // This member function returns one of the infinitely many
    // points on the plane that could serve as the origin of a 2D
    // coordinate system.
    template <class Type>
    brick::numeric::Vector3D<Type> const&
    Plane3D<Type>::
    getOrigin() const
    {
      return m_origin;
    }


    template <class Type>
    Type
    Plane3D<Type>::
    findDistance(brick::numeric::Vector3D<Type> const& point) const
    {
      brick::numeric::Vector3D<Type> offset = point - m_origin;
      offset -= (brick::numeric::dot<Type>(offset, m_directionVector0)
                 * m_directionVector0);
      offset -= (brick::numeric::dot<Type>(offset, m_directionVector1)
                 * m_directionVector1);
      return brick::numeric::magnitude<Type>(offset);
    }


    /* ======= Private member functions. ======= */

    template <class Type>
    template <class Iterator>
    void
    Plane3D<Type>::
    estimateFromSequence(Iterator beginIterator,
                         Iterator endIterator)
    {
      // Get mean point.
      brick::numeric::Vector3D<Type> meanPoint;
      Iterator targetIterator = beginIterator;
      size_t count = 0;
      while(targetIterator != endIterator) {
        meanPoint += *targetIterator;
        ++count;
        ++targetIterator;
      }
      meanPoint /= static_cast<double>(count);

      // Get 3x3 covariance matrix.
      brick::numeric::Array2D<double> covarianceMatrix(3, 3);
      covarianceMatrix = 0.0;
      targetIterator = beginIterator;
      while(targetIterator != endIterator) {
        const double xValue = targetIterator->x() - meanPoint.x();
        const double yValue = targetIterator->y() - meanPoint.y();
        const double zValue = targetIterator->z() - meanPoint.z();
        covarianceMatrix(0, 0) += xValue * xValue;
        covarianceMatrix(0, 1) += xValue * yValue;
        covarianceMatrix(0, 2) += xValue * zValue;
        covarianceMatrix(1, 1) += yValue * yValue;
        covarianceMatrix(1, 2) += yValue * zValue;
        covarianceMatrix(2, 2) += zValue * zValue;
        ++targetIterator;
      }
      covarianceMatrix(1, 0) = covarianceMatrix(0, 1);
      covarianceMatrix(2, 0) = covarianceMatrix(0, 2);
      covarianceMatrix(2, 1) = covarianceMatrix(1, 2);
      covarianceMatrix /= static_cast<double>(count);

      // Solve for best fit plane.
      brick::numeric::Array1D<double> eigenvalues;
      brick::numeric::Array2D<double> eigenvectors;
      brick::linearAlgebra::eigenvectorsSymmetric(
        covarianceMatrix, eigenvalues, eigenvectors);
      brick::numeric::Vector3D<Type> direction0(
        eigenvectors(0, 0), eigenvectors(1, 0), eigenvectors(2, 0));
      brick::numeric::Vector3D<Type> direction1(
        eigenvectors(0, 1), eigenvectors(1, 1), eigenvectors(2, 1));
      direction0 /= brick::numeric::magnitude<Type>(direction0);
      direction1 /= brick::numeric::magnitude<Type>(direction1);

      m_origin = meanPoint;
      m_directionVector0 = direction0;
      m_directionVector1 = direction1;
    }


    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Plane3D<Type>& plane)
    {
      stream << "Plane3D{ "
             << plane.getOrigin() << ", "
             << plane.getDirectionVector0() << ", "
             << plane.getDirectionVector1() << " }";
      return stream;
    }


  } // namespace utilities

} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_PLANE3D_IMPL_HH */
