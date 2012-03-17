/**
***************************************************************************
* @file brick/computerVision/calibrationTools_impl.hh
*
* Header file declaring utility routines for use in calibrating sensors.
*
* Copyright (C) 2009-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CALIBRATIONTOOLS_IMPL_HH
#define BRICK_COMPUTERVISION_CALIBRATIONTOOLS_IMPL_HH

// This file is included by calibrationTools.hh, and should not be
// directly included by user code, so no need to include
// calibrationTools.hh here.
// 
// #include <brick/computerVision/calibrationTools.hh>

#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/rotations.hh>
#include <brick/numeric/transform3D.hh>
#include <brick/numeric/transform3DTo2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/optimization/gradientFunctionLM.hh>
#include <brick/optimization/optimizerLM.hh>

namespace brick {

  namespace computerVision {

    namespace privateCode {

      /**
       ** This class implements the the sum-of-squares error function
       ** that will be minimized by estimateCameraIntrinsics.
       **/
      template<class Intrinsics>
      class CameraIntrinsicsObjectiveFunction
        : public std::unary_function<
            typename Intrinsics::ParameterVectorType,
            brick::numeric::Array1D<typename Intrinsics::FloatType> >
      {
      public:

        /// Convenient typedef for keeping track of how we represent
        /// real numbers.
        typedef typename Intrinsics::FloatType FloatType;
        
        /**
         * Constructor.
         * 
         * @param intrinsics This argument gives initial values for
         * the state of *this.  It does not have to have all if its
         * intrinsic parameters set correctly, but it should
         * accurately reflect the image size (numPixelsX, numPixelsY).
         * 
         * @param points3DBegin This argument is the beginning of a
         * sequence of Vector3D instances expressed in camera
         * coordinates.
         * 
         * @param points3DEnd This argument This argument marks the
         * end of the sequence of Vector3D instances expressed in
         * camera coordinates.
         * 
         * @param points2DBegin This argument is the beginning of a
         * sequence of Vector2D instances in image coordinates,
         * corresponding to the 3D sequence in earlier constructor
         * arguments.
         */
        template <class Iter3D, class Iter2D>
        CameraIntrinsicsObjectiveFunction(Intrinsics intrinsics,
                                          Iter3D points3DBegin,
                                          Iter3D points3DEnd,
                                          Iter2D points2DBegin)
          : m_cameraCoords(points3DBegin, points3DEnd),
            m_imageCoords(points2DBegin, points2DBegin + m_cameraCoords.size()),
            m_intrinsics(intrinsics),
            m_numPixelsX(intrinsics.getNumPixelsX()),
            m_numPixelsY(intrinsics.getNumPixelsY())
          {}


        /**
         * Constructor for use by subclasses.
         * 
         * @param intrinsics This argument gives initial values for
         * the state of *this.  It does not have to have all if its
         * intrinsic parameters set correctly, but it should
         * accurately reflect the image size (numPixelsX, numPixelsY).
         */
        CameraIntrinsicsObjectiveFunction(Intrinsics const& intrinsics)
          : m_cameraCoords(),
            m_imageCoords(),
            m_intrinsics(intrinsics),
            m_numPixelsX(intrinsics.getNumPixelsX()),
            m_numPixelsY(intrinsics.getNumPixelsY())
          {}
        
        
        /**
         * Destructor.
         */
        virtual
        ~CameraIntrinsicsObjectiveFunction() {}


        /**
         * Application operator returns a sequence of error terms that
         * should be squared and added together.
         *
         * @param freeParameters This argument specifies the minimal
         * subset of camera extrinsics & intrinsics to be evaluated.
         * 3D points will be projected based on these parameters and
         * the image-space residuals will be returned.
         *
         * @return The return value is an array of image-space
         * residuals.
         */
        virtual numeric::Array1D<FloatType>
        operator()(typename Intrinsics::ParameterVectorType const&
                   freeParameters) {
          this->setParameters(freeParameters);
          numeric::Array1D<FloatType> errorTerms(m_cameraCoords.size() * 2);
          for(unsigned int ii = 0; ii < m_cameraCoords.size(); ++ii) {
            numeric::Vector2D<FloatType> projectedPoint = m_intrinsics.project(
              m_cameraCoords[ii]);
            numeric::Vector2D<FloatType> residualVector =
              m_imageCoords[ii] - projectedPoint;
            errorTerms[2 * ii] = residualVector.x();
            errorTerms[(2 * ii) + 1] = residualVector.y();
          }
          return errorTerms;
        }


        virtual Intrinsics
        getIntrinsics() {return m_intrinsics;}


        /** 
         * Returns the number of 3D points in the sequence provided to
         * the constructor..
         * 
         * @return The return value the number of points.
         */
        virtual unsigned int
        getNumberOfPoints() {return m_cameraCoords.size();}
        

        /** 
         * Updates the internal state to reflect the specified
         * parameter vector.  This function differs from
         * this->operator()() in that it stops before computing any
         * image-space residuals.
         * 
         * @param parameters This argument exactly corresponds to the
         * "freeParameters" argument of this->operator()().
         */
        virtual void
        setParameters(typename Intrinsics::ParameterVectorType const&
                      parameters) {
          m_intrinsics.setFreeParameters(parameters);
          this->estimateDependentParameters(
            m_intrinsics, m_numPixelsX, m_numPixelsY,
            m_cameraCoords.begin(), m_cameraCoords.end(),
            m_imageCoords.begin());
        }

      protected:

        // Given intrinsics with known lens distortion and point
        // correspondences between camera coords and image coords,
        // solve for the remaining (pinhole) intrinsics using
        // estimateCameraIntrinsicsPinhole().
        template<class Iter3D, class Iter2D>
        void
        estimateDependentParameters(Intrinsics& intrinsics,
                                    unsigned int numPixelsX, unsigned int numPixelsY,
                                    Iter3D points3DBegin, Iter3D points3DEnd,
                                    Iter2D points2DBegin) {
          std::vector< numeric::Vector3D<FloatType> > undistortedPoints(
            points3DEnd - points3DBegin);
          for(unsigned int ii = 0; ii < undistortedPoints.size(); ++ii) {
            undistortedPoints[ii] = intrinsics.projectThroughDistortion(
              *points3DBegin);
            ++points3DBegin;
          }
          CameraIntrinsicsPinhole<FloatType> dependentIntrinsics =
            estimateCameraIntrinsicsPinhole<FloatType>(
              numPixelsX, numPixelsY, undistortedPoints.begin(),
              undistortedPoints.end(), points2DBegin);
          intrinsics.setDependentParameters(
            dependentIntrinsics.getNumPixelsX(),
            dependentIntrinsics.getNumPixelsY(),
            (dependentIntrinsics.getFocalLength()
             / dependentIntrinsics.getPixelSizeX()),
            (dependentIntrinsics.getFocalLength()
             / dependentIntrinsics.getPixelSizeY()),
            dependentIntrinsics.getCenterU(),
            dependentIntrinsics.getCenterV());
        }

        std::vector< brick::numeric::Vector3D<FloatType> > m_cameraCoords;
        std::vector< brick::numeric::Vector2D<FloatType> > m_imageCoords;
        Intrinsics m_intrinsics;
        unsigned int m_numPixelsX;
        unsigned int m_numPixelsY;
      };


      /**
       ** This class implements the the sum-of-squares error function
       ** that will be minimized by estimateCameraParameters.
       **/
      template<class Intrinsics>
      class CameraParametersObjectiveFunction
        : public CameraIntrinsicsObjectiveFunction<Intrinsics>
      {
      public:

        /// Convenient typedef for keeping track of how we represent
        /// real numbers.
        typedef typename Intrinsics::FloatType FloatType;

        /**
         * Constructor.
         * 
         * @param intrinsics This argument gives initial values for
         * the state of *this.  It does not have to have all if its
         * intrinsic parameters set correctly, but it should
         * accurately reflect the image size (numPixelsX, numPixelsY).
         * 
         * @param points3DBegin This argument is the beginning of a
         * sequence of Vector3D instances expressed in world
         * coordinates (not camera coordinates).
         * 
         * @param points3DEnd This argument This argument marks the
         * end of the sequence of Vector3D instances expressed in
         * world coordinates.
         * 
         * @param points2DBegin This argument is the beginning of a
         * sequence of Vector2D instances in image coordinates,
         * corresponding to the 3D sequence in earlier constructor
         * arguments.
         */
        template <class Iter3D, class Iter2D>
        CameraParametersObjectiveFunction(
          Intrinsics const& intrinsics,
          Iter3D points3DBegin, Iter3D points3DEnd,
          Iter2D points2DBegin)
          : CameraIntrinsicsObjectiveFunction<Intrinsics>(intrinsics),
            m_intrinsicParameters(intrinsics.getNominalFreeParameters()),
            m_worldCoords(points3DBegin, points3DEnd),
            m_cameraTworld()
          {
            this->m_cameraCoords.resize(m_worldCoords.size());
            this->m_imageCoords.resize(m_worldCoords.size());
            std::copy(points2DBegin, points2DBegin + m_worldCoords.size(),
                      this->m_imageCoords.begin());
          }
        
        
        /**
         * Destructor.
         */
        virtual
        ~CameraParametersObjectiveFunction() {};


        virtual numeric::Transform3D<FloatType>
        getPoseCameraTworld() {return m_cameraTworld;}

        
        virtual void
        setParameters(typename Intrinsics::ParameterVectorType const&
                      parameters) {
          // Update camera coords to reflect extrinsics parameters.
          numeric::Quaternion<FloatType> quaternion(
            parameters[0], parameters[1], parameters[2], parameters[3]);
          quaternion.normalize();
          m_cameraTworld = numeric::quaternionToTransform3D(quaternion);
          m_cameraTworld.setValue(0, 3, parameters[4]);
          m_cameraTworld.setValue(1, 3, parameters[5]);
          m_cameraTworld.setValue(2, 3, parameters[6]);
          std::transform(m_worldCoords.begin(), m_worldCoords.end(),
                         this->m_cameraCoords.begin(),
                         m_cameraTworld.getFunctor());

          // Parent class doesn't know about extrinsics, so we have to
          // repackage intrinsic parameters before we can dispatch to
          // parent's version of setParameters.
          std::copy(parameters.begin() + 7, parameters.end(),
                    m_intrinsicParameters.begin());
          CameraIntrinsicsObjectiveFunction<Intrinsics>::setParameters(
            m_intrinsicParameters);
        }

      private:

        // Members solely to avoid reallocating resources on each
        // function call.  Otherwise non-essential.
        typename Intrinsics::ParameterVectorType m_intrinsicParameters;

        // New members need for this subclass.
        std::vector< numeric::Vector3D<FloatType> > m_worldCoords;
        numeric::Transform3D<FloatType> m_cameraTworld;
      };


      // Useful for figuring out whether points are in of the camera
      // or behind it.
      template<class FloatType, class Iter3D>
      FloatType
      computeAverageTransformedZValue(
        numeric::Transform3D<FloatType> const& cameraTworld,
        Iter3D points3DBegin, Iter3D points3DEnd)
      {
        FloatType zValue = 0.0;
        unsigned int count = 0;
        while(points3DBegin != points3DEnd) {
          Vector3D<FloatType> transformedPoint =
            cameraTworld * (*points3DBegin);
          zValue += transformedPoint.z();
          ++count;
          ++points3DBegin;
        }
        if(count != 0) {
          zValue /= count;
        }
        return zValue;
      }
      

#if 0
      // xxx save this because it'll be useful for adapting to
      // estimateTransform3D and estimateTransform2D.
      template <class FloatType, class Iter3D, class Iter2D>
      numeric::Transform3DTo2D<FloatType>
      estimateTransform3DTo2D(Iter3D points3DBegin, Iter3D points3DEnd,
                              Iter2D points2DBegin)
      {
        // If points3D are (x_i, y_i, z_i) and points2D are (u_i,
        // v_i), this routine returns a Transform3DTo2D H so that for
        // each point in points3D and corresponding point in points2D:
        // 
        //   u_i ~= ((h0*x_i + h1*y_i + h2*z_i + h3)
        //           / (h8*x_i + h9*y_i + h10*z_i + 1.0))
        // 
        //   v_i ~= ((h4*x_i + h5*y_i + h6*z_i + h7)
        //           / (h8*x_i + h9*y_i + h10*z_i + 1.0))
        // 
        // where
        // 
        //        [[h0, h1, h2,   h3],
        //   H ==  [h4, h5, h6,   h7],
        //         [h8, h9, h10, 1.0]]
        // 
        // A linear solution is used, however the equations are scaled by
        // the denominators of the above equations, so the solution is _not_
        // optimal in the least squares sense.  If iter == 0, the solution is
        // refined iter times.  At each iteration, the equations are rescaled
        // based on the previous estimate of H.

        // Depending on the type of iterator, this might be O(N), but the
        // rest of the function will dominate regardless.
        unsigned int numberOfInputPoints = points3DEnd - points3DBegin;

        numeric::Array2D<FloatType> AMatrix(2 * numberOfInputPoints, 11);
        numeric::Array1D<FloatType> bVector(2 * numberOfInputPoints);
        unsigned int rowIndex = 0;
        while(points3DBegin != points3DEnd) {
          AMatrix(rowIndex, 0) = points3DBegin->x();
          AMatrix(rowIndex, 1) = points3DBegin->y();
          AMatrix(rowIndex, 2) = points3DBegin->z();
          AMatrix(rowIndex, 3) = 1.0;
          AMatrix(rowIndex, 4) = 0.0;
          AMatrix(rowIndex, 5) = 0.0;
          AMatrix(rowIndex, 6) = 0.0;
          AMatrix(rowIndex, 7) = 0.0;
          AMatrix(rowIndex, 8) = -(points3DBegin->x() * points2DBegin->x());
          AMatrix(rowIndex, 9) = -(points3DBegin->y() * points2DBegin->x());
          AMatrix(rowIndex, 10) = -(points3DBegin->z()
                                    * points2DBegin->x());
          bVector[rowIndex] = points2DBegin->x();

          AMatrix(rowIndex + 1, 0) = 0.0;
          AMatrix(rowIndex + 1, 1) = 0.0;
          AMatrix(rowIndex + 1, 2) = 0.0;
          AMatrix(rowIndex + 1, 3) = 0.0;
          AMatrix(rowIndex + 1, 4) = points3DBegin->x();
          AMatrix(rowIndex + 1, 5) = points3DBegin->y();
          AMatrix(rowIndex + 1, 6) = points3DBegin->z();
          AMatrix(rowIndex + 1, 7) = 1.0;
          AMatrix(rowIndex + 1, 8) = -(points3DBegin->x() * points2DBegin->y());
          AMatrix(rowIndex + 1, 9) = -(points3DBegin->y() * points2DBegin->y());
          AMatrix(rowIndex + 1, 10) = -(points3DBegin->z()
                                        * points2DBegin->y());
          bVector[rowIndex + 1] = points2DBegin->y();
          rowIndex += 2;
          ++points2DBegin;
          ++points3DBegin;
        }

        // Now we can solve for our unknown parameters using the
        // Moore-Penrose equations.
        numeric::Array2D<FloatType> ATranspose = AMatrix.transpose();
        numeric::Array2D<FloatType> ATA = numeric::matrixMultiply<FloatType>(
          ATranspose, AMatrix);
        numeric::Array1D<FloatType> ATb = numeric::matrixMultiply<FloatType>(
          ATranspose, bVector);
        linearAlgebra::linearSolveInPlace(ATA, ATb);

        // Result is returned in ATb.
        numeric::Transform3DTo2D<FloatType> result(
          ATb[0], ATb[1],  ATb[2],  ATb[3],
          ATb[4], ATb[5],  ATb[6],  ATb[7],
          ATb[8], ATb[9],  ATb[10],  1.0);
        return result;
      }

#endif /* #if 0 */

      
      template <class FloatType, class Iter3D, class Iter2D>
      numeric::Transform3DTo2D<FloatType>
      estimateTransform3DTo2D(Iter3D points3DBegin, Iter3D points3DEnd,
                              Iter2D points2DBegin)
      {
        // If points3D are (x_i, y_i, z_i) and points2D are (u_i,
        // v_i), this routine returns a Transform3DTo2D H so that for
        // each point in points3D and corresponding point in points2D:
        // 
        //   u_i ~= ((h0*x_i + h1*y_i + h2*z_i + h3)
        //           / (h8*x_i + h9*y_i + h10*z_i + h11))
        // 
        //   v_i ~= ((h4*x_i + h5*y_i + h6*z_i + h7)
        //           / (h8*x_i + h9*y_i + h10*z_i + h11))
        // 
        // where
        // 
        //        [[h0, h1,  h2,  h3],
        //   H ==  [h4, h5,  h6,  h7],
        //         [h8, h9, h10, h11]]
        // 
        // A linear solution is used, however the equations are scaled by
        // the denominators of the above equations, so the solution is _not_
        // optimal in the least squares sense.  If iter == 0, the solution is
        // refined iter times.  At each iteration, the equations are rescaled
        // based on the previous estimate of H.

        // Depending on the type of iterator, this might be O(N), but the
        // rest of the function will dominate regardless.
        unsigned int numberOfInputPoints = points3DEnd - points3DBegin;

        numeric::Array2D<FloatType> AMatrix(2 * numberOfInputPoints, 12);
        unsigned int rowIndex = 0;
        while(points3DBegin != points3DEnd) {
          AMatrix(rowIndex, 0) = points3DBegin->x();
          AMatrix(rowIndex, 1) = points3DBegin->y();
          AMatrix(rowIndex, 2) = points3DBegin->z();
          AMatrix(rowIndex, 3) = 1.0;
          AMatrix(rowIndex, 4) = 0.0;
          AMatrix(rowIndex, 5) = 0.0;
          AMatrix(rowIndex, 6) = 0.0;
          AMatrix(rowIndex, 7) = 0.0;
          AMatrix(rowIndex, 8) = -(points3DBegin->x() * points2DBegin->x());
          AMatrix(rowIndex, 9) = -(points3DBegin->y() * points2DBegin->x());
          AMatrix(rowIndex, 10) = -(points3DBegin->z()
                                    * points2DBegin->x());
          AMatrix(rowIndex, 11) = -points2DBegin->x();

          AMatrix(rowIndex + 1, 0) = 0.0;
          AMatrix(rowIndex + 1, 1) = 0.0;
          AMatrix(rowIndex + 1, 2) = 0.0;
          AMatrix(rowIndex + 1, 3) = 0.0;
          AMatrix(rowIndex + 1, 4) = points3DBegin->x();
          AMatrix(rowIndex + 1, 5) = points3DBegin->y();
          AMatrix(rowIndex + 1, 6) = points3DBegin->z();
          AMatrix(rowIndex + 1, 7) = 1.0;
          AMatrix(rowIndex + 1, 8) = -(points3DBegin->x() * points2DBegin->y());
          AMatrix(rowIndex + 1, 9) = -(points3DBegin->y() * points2DBegin->y());
          AMatrix(rowIndex + 1, 10) = -(points3DBegin->z()
                                        * points2DBegin->y());
          AMatrix(rowIndex + 1, 11) = -points2DBegin->y();
          rowIndex += 2;
          ++points2DBegin;
          ++points3DBegin;
        }

        // Now we can solve for our unknown parameters using the
        // smallest eigenvalue of AMatrix^T * AMatrix.
        numeric::Array2D<FloatType> ATranspose = AMatrix.transpose();
        numeric::Array2D<FloatType> ATA = numeric::matrixMultiply<FloatType>(
          ATranspose, AMatrix);
        numeric::Array1D<FloatType> eigenvalues;
        numeric::Array2D<FloatType> eigenvectors;
        linearAlgebra::eigenvectorsSymmetric(ATA, eigenvalues, eigenvectors);

        // Unpack the result.
        unsigned int lastColumn = eigenvectors.columns() - 1;
        numeric::Transform3DTo2D<FloatType> result(
          eigenvectors(0, lastColumn), eigenvectors(1, lastColumn),
          eigenvectors(2, lastColumn), eigenvectors(3, lastColumn),
          eigenvectors(4, lastColumn), eigenvectors(5, lastColumn),
          eigenvectors(6, lastColumn), eigenvectors(7, lastColumn),
          eigenvectors(8, lastColumn), eigenvectors(9, lastColumn),
          eigenvectors(10, lastColumn), eigenvectors(11, lastColumn));
        return result;
      }
      
    } // namespace privateCode

    
    // This function estimates camera intrinsic parameters for
    // "complicated" types of camera intrinsics that require nonlinear
    // optimization.
    template <class Intrinsics, class Iter3D, class Iter2D>
    Intrinsics
    estimateCameraIntrinsics(unsigned int numPixelsX, unsigned int numPixelsY,
                             Iter3D points3DBegin,
                             Iter3D points3DEnd,
                             Iter2D points2DBegin,
                             int verbosity)
    {
      typedef privateCode::CameraIntrinsicsObjectiveFunction<Intrinsics>
        ObjectiveFunction;
      typedef brick::optimization::GradientFunctionLM<ObjectiveFunction>
        GradientFunctionLM;
      
      Intrinsics intrinsics;
      intrinsics.setNumPixelsX(numPixelsX);
      intrinsics.setNumPixelsY(numPixelsY);
      typename Intrinsics::ParameterVectorType freeParameters =
        intrinsics.getNominalFreeParameters();
      ObjectiveFunction objectiveFunction(
        intrinsics, points3DBegin, points3DEnd, points2DBegin);

      GradientFunctionLM gradientFunction(objectiveFunction);
      OptimizerLM<GradientFunctionLM> optimizer(gradientFunction);
      optimizer.setMaxIterations(100);
      optimizer.setMaxLambda(1.0E15);
      optimizer.setMinDrop(1.0E-6);
      optimizer.setVerbosity(verbosity);
      optimizer.setStartPoint(freeParameters);
      freeParameters = optimizer.optimum();

      objectiveFunction.setParameters(freeParameters);
      return objectiveFunction.getIntrinsics();
    }


#if 0  /* Keep this around temporarily because it's so pretty. */
    // This function estimates pinhole camera intrinsic parameters
    // based on corresponding points in 2D image coordinates and 3D
    // camera coordinates.
    template <FloatType, class Iter3D, class Iter2D>
    CameraIntrinsicsPinhole<FloatType>
    estimateCameraIntrinsicsPinhole(
      unsigned int numPixelsX, unsigned int numPixelsY,
      Iter3D points3DBegin, Iter3D points3DEnd,
      Iter2D points2DBegin)
    {
      // Depending on the type of iterator, this might be O(N), but the
      // rest of the function will dominate regardless.
      unsigned int numberOfInputPoints = points3DEnd - points3DBegin;

      // For each 2D point, w_i = [u_i, v_i, 1]^T, expressed in 2D
      // homogeneous coordinates, and corresponding 3D point, d_i =
      // [x_i, y_i, z_i, 1]^T, expressed in 3D homogeneous coordinates,
      // the projection equation is
      //
      //   alpha * w_i = P * d_i
      // 
      // where alpha is an arbitrary scale factor, and P is the pinhole
      // projection matrix.
      // 
      //       | k_x, 0.0, u_0, 0.0 |
      //   P = | 0.0, k_y, v_0, 0.0 |
      //       | 0.0, 0.0, 1.0, 0.0 |
      //
      // This implies that w_i and (P * d_i) are parallel,
      // which implies that their cross product has zero magnitude:
      //
      //   cross(w_i, P * d_i) = [0.0, 0.0, 0.0]^T
      //
      // This cross product can be computed as follows (the cross()
      // and skewSymmetric() functions are declared in
      // brick/numeric/utilities.hh):
      //
      //   cross(w_i, P * d_i) = skewSymmetric(w_i) * (P * d_i)
      //
      // We can rearrange the matrix*vector product P * d_i to get
      // Equation 1:
      //
      //   cross(w_i, P * d_i)
      //     = skewSymmetric(w_i) * equivalentMatrix(d_i) * vec(P)
      //     = [0.0, 0.0, 0.0]^T
      //    
      // We get one of these sets of linear equations in the elements of
      // P for each pair of input points.  We can assemble the equations
      // into a linear system and solve for vec(P).
      // 
      //   A * vec(P) = [0.0, 0.0, ..., 0.0]
      //
      // In practice, though, we see that P is very sparse, and has one
      // element identically equal to one.  The first of these two
      // observations means that we must drop many of the columns of
      // equivalentMatrix(d_i).  Rather than use the prefab function
      // "equivalentMatrix()" and ignore or remove unused columns, we
      // hand-roll our own sparse equivalentMatrix below.
      //
      // The second observation (that the [2, 3] element of P is equal
      // to 1.0) lets us rearrange the matrix equation above.  Here is
      // the P * d_i product (the equivalentMatrix(d_i) * vec(P)
      // product), written out by individual elements, retaining only
      // the nonzero elements of vec(P):
      //
      //   | x_i, 0.0, z_i, 0.0 |   | k_x |   | 0.0 |
      //   | 0.0, y_i, 0.0, z_i | * | k_y | + | 0.0 | 
      //   | 0.0, 0.0, 0.0, 0.0 |   | u_0 |   | z_i |
      //                            | v_0 |
      //
      // The added term, [0.0, 0.0, z_i]^T, corresponds to the [2, 3]
      // element of P.
      // 
      // Here is skewSymmetric(w_i) written out:
      //
      //                        |  0.0, -1.0,  v_i |
      //   skewSymmetric(w_i) = |  1.0,  0.0, -u_i |
      //                        | -v_i,  u_i,  0.0 |
      //
      // Substituting the previous two equations into Eq. 1, we have:
      //
      //   |      0.0,    -y_i,      0.0,   -z_i |   | k_x |   |  v_i * z_i |
      //   |      x_i,     0.0,      z_i,    0.0 | * | k_y | + | -u_i * z_i |
      //   | -v_i*x_i, u_i*y_i, -v_i*z_i, u_i*z_i|   | u_0 |   |     0.0    |
      //                                             | v_0 |
      //
      //          | 0.0 |
      //        = | 0.0 |
      //          | 0.0 |
      //
      // Subtracting from both sides the vector that does not depend on
      // our unknown parameters gives:
      //
      //   |      0.0,    -y_i,      0.0,   -z_i |   | k_x |   | -v_i * z_i |
      //   |      x_i,     0.0,      z_i,    0.0 | * | k_y | = |  u_i * z_i |
      //   | -v_i*x_i, u_i*y_i, -v_i*z_i, u_i*z_i|   | u_0 |   |     0.0    |
      //                                             | v_0 |
      // These are the equations implemented below.
      numeric::Array2D<FloatType> AMatrix(3 * numberOfInputPoints, 4);
      numeric::Array1D<FloatType> bVector(3 * numberOfInputPoints);
      unsigned int rowIndex = 0;
      while(points3DBegin != points3DEnd) {
        AMatrix(rowIndex, 0) = 0.0;
        AMatrix(rowIndex, 1) = -(points3DBegin->y());
        AMatrix(rowIndex, 2) = 0.0;
        AMatrix(rowIndex, 3) = -(points3DBegin->z());
        AMatrix(rowIndex + 1, 0) = points3DBegin->x();
        AMatrix(rowIndex + 1, 1) = 0.0;
        AMatrix(rowIndex + 1, 2) = points3DBegin->z();
        AMatrix(rowIndex + 1, 3) = 0.0;
        AMatrix(rowIndex + 2, 0) = -(points2DBegin->y()) * points3DBegin->x();
        AMatrix(rowIndex + 2, 1) = points2DBegin->x() * points3DBegin->y();
        AMatrix(rowIndex + 2, 2) = -(points2DBegin->y()) * points3DBegin->z();
        AMatrix(rowIndex + 2, 3) = points2DBegin->x() * points3DBegin->z();
        bVector(rowIndex) = AMatrix(rowIndex + 2, 2);
        bVector(rowIndex + 1) = AMatrix(rowIndex + 2, 3);
        bVector(rowIndex + 2) = 0.0;

        rowIndex += 3;
        ++points2DBegin;
        ++points3DBegin;
      }

      // Now we can solve for our unknown parameters using the
      // Moore-Penrose equations.
      numeric::Array2D<FloatType> ATranspose = AMatrix.transpose();
      numeric::Array2D<FloatType> ATA = numeric::matrixMultiply<FloatType>(
        ATranspose, AMatrix);
      numeric::Array1D<FloatType> ATb = numeric::matrixMultiply<FloatType>(
        ATranspose, bVector);
      linearAlgebra::linearSolveInPlace(ATA, ATb);

      // Result is returned in ATb.
      FloatType k_x = ATb[0];
      FloatType k_y = ATb[1];
      FloatType u_0 = ATb[2];
      FloatType v_0 = ATb[3];

      return CameraIntrinsicsPinhole(
        numPixelsX, numPixelsY, 1, 1.0 / k_x, 1.0 / k_y, u_0, v_0);
    }
#endif /* #if 0 */


    // This function estimates pinhole camera intrinsic parameters
    // based on corresponding points in 2D image coordinates and 3D
    // camera coordinates.
    template <class FloatType, class Iter3D, class Iter2D>
    CameraIntrinsicsPinhole<FloatType>
    estimateCameraIntrinsicsPinhole(unsigned int numPixelsX,
                                    unsigned int numPixelsY,
                                    Iter3D points3DBegin,
                                    Iter3D points3DEnd,
                                    Iter2D points2DBegin)
    {
      // Depending on the type of iterator, this might be O(N), but the
      // rest of the function will dominate regardless.
      unsigned int numberOfInputPoints = points3DEnd - points3DBegin;

      // For each 2D point, w_i = [u_i, v_i, 1]^T, expressed in 2D
      // homogeneous coordinates, and corresponding 3D point, d_i =
      // [x_i, y_i, z_i, 1]^T, expressed in 3D homogeneous coordinates,
      // the projection equation is
      //
      //   alpha * w_i = P * d_i
      // 
      // where alpha is an arbitrary scale factor, and P is the pinhole
      // projection matrix.
      // 
      //       | k_x, 0.0, u_0, 0.0 |
      //   P = | 0.0, k_y, v_0, 0.0 |
      //       | 0.0, 0.0, 1.0, 0.0 |
      //
      // Note that, by this definition, alpha is always equal to z_i:
      // points that lie on the same ray in camera coordinates project
      // to the same image point.  Dividing out the alpha, we have
      // 
      //   w_i = P * [x_i / z_i, y_i / z_i, 1, 1]^T
      //
      // We're going to pick elements of P to minimize (in the least
      // squares sense) the differences between the right and left
      // sides of this equation for all the input points.  Rescaling
      // by z_i has the added advantage of normalizing the
      // contributions of the input points so that we're essentially
      // minimizing an image-space residual.
      //
      // Expanding the above equation, we have
      //
      //   |u_i|   | k_x, 0.0, u_0, 0.0 |   |x_i / z_i|
      //   |v_i| = | 0.0, k_y, v_0, 0.0 | * |y_i / z_i| 
      //   | 1 |   | 0.0, 0.0, 1.0, 0.0 |   |    1    |
      //                                    |    1    |
      // The bottom row of this equation is an identity (1 = 1), so we
      // discard it and rearrange the top two rows to isolate the unknowns.
      // 
      //   |u_i| = | k_x * x_i / z_i| + |u_0|
      //   |v_i|   | k_y * y_i / z_i|   |v_0|
      //
      //   |x_i / z_i,     0    , 1, 0|   |k_x|   [u_0|
      //   |    0    , y_i / z_i, 0, 1| * |k_y| = |v_0|
      //                                  |u_0|
      //                                  |v_0|
      // 
      // Combining these equations for all of the input points, we
      // solve for k_x, k_y, u_0, and v_0 below.
      numeric::Array2D<FloatType> AMatrix(2 * numberOfInputPoints, 4);
      numeric::Array1D<FloatType> bVector(2 * numberOfInputPoints);
      unsigned int rowIndex = 0;
      while(points3DBegin != points3DEnd) {
        AMatrix(rowIndex, 0) = points3DBegin->x() / points3DBegin->z();
        AMatrix(rowIndex, 1) = 0.0;
        AMatrix(rowIndex, 2) = 1.0;
        AMatrix(rowIndex, 3) = 0.0;
        AMatrix(rowIndex + 1, 0) = 0.0;
        AMatrix(rowIndex + 1, 1) = points3DBegin->y() / points3DBegin->z();
        AMatrix(rowIndex + 1, 2) = 0.0;
        AMatrix(rowIndex + 1, 3) = 1.0;
        bVector(rowIndex) = points2DBegin->x();
        bVector(rowIndex + 1) = points2DBegin->y();
        rowIndex += 2;
        ++points2DBegin;
        ++points3DBegin;
      }

      // Now we can solve for our unknown parameters using the
      // Moore-Penrose equations.
      numeric::Array2D<FloatType> ATranspose = AMatrix.transpose();
      numeric::Array2D<FloatType> ATA = numeric::matrixMultiply<FloatType>(
        ATranspose, AMatrix);
      numeric::Array1D<FloatType> ATb = numeric::matrixMultiply<FloatType>(
        ATranspose, bVector);
      linearAlgebra::linearSolveInPlace(ATA, ATb);

      // Result is returned in ATb.
      FloatType k_x = ATb[0];
      FloatType k_y = ATb[1];
      FloatType u_0 = ATb[2];
      FloatType v_0 = ATb[3];

      return CameraIntrinsicsPinhole<FloatType>(
        numPixelsX, numPixelsY, 1, 1.0 / k_x, 1.0 / k_y, u_0, v_0);
    }
    

    template <class Intrinsics, class Iter3D, class Iter2D>
    void
    estimateCameraParameters(
      Intrinsics& intrinsics,
      numeric::Transform3D<typename Intrinsics::FloatType>& cameraTworld,
      unsigned int numPixelsX, unsigned int numPixelsY,
      Iter3D points3DBegin, Iter3D points3DEnd,
      Iter2D points2DBegin,
      int verbosity)
    {
      typedef privateCode::CameraParametersObjectiveFunction<Intrinsics>
        ObjectiveFunction;
      typedef brick::optimization::GradientFunctionLM<ObjectiveFunction>
        GradientFunctionLM;
      typedef typename Intrinsics::FloatType FloatType;

      // Set output values that we know already.
      intrinsics.setNumPixelsX(numPixelsX);
      intrinsics.setNumPixelsY(numPixelsY);

      // Get rough initial guess at pinhole intrinsics and extrinsics.
      // We'll use this to refine the starting point for our
      // optimization.
      numeric::Transform3D<FloatType> camTworldInitial;
      CameraIntrinsicsPinhole<FloatType> initialIntrinsics;
      estimateCameraParametersPinhole(
        initialIntrinsics, camTworldInitial,
        numPixelsX, numPixelsY, points3DBegin, points3DEnd, points2DBegin);

      // Get a nominal starting point for the lens distortion parameters.
      typename Intrinsics::ParameterVectorType intrinsicFreeParameters =
        intrinsics.getNominalFreeParameters();

      // Combine initial pinhole parameters, extrinsics, and lens
      // distortion into our the starting parameter vector for our
      // optimization.
      typename Intrinsics::ParameterVectorType allFreeParameters(
        intrinsicFreeParameters.size() + 7);
      numeric::Quaternion<FloatType> camQworld =
        numeric::transform3DToQuaternion(camTworldInitial);
      allFreeParameters[0] = camQworld.s();
      allFreeParameters[1] = camQworld.i();
      allFreeParameters[2] = camQworld.j();
      allFreeParameters[3] = camQworld.k();
      allFreeParameters[4] = camTworldInitial(0, 3);
      allFreeParameters[5] = camTworldInitial(1, 3);
      allFreeParameters[6] = camTworldInitial(2, 3);
      std::copy(intrinsicFreeParameters.begin(), intrinsicFreeParameters.end(),
                allFreeParameters.begin() + 7);

      // Run the optimization.
      ObjectiveFunction objectiveFunction(
        intrinsics, points3DBegin, points3DEnd, points2DBegin);
      GradientFunctionLM gradientFunction(objectiveFunction);
      OptimizerLM<GradientFunctionLM> optimizerLM(gradientFunction);
      optimizerLM.setMaxIterations(100);
      optimizerLM.setMaxLambda(1.0E15);
      optimizerLM.setMinDrop(1.0E-6);
      optimizerLM.setVerbosity(verbosity);
      optimizerLM.setStartPoint(allFreeParameters);
      allFreeParameters = optimizerLM.optimum();

      // Communicate the result back to the calling context.
      objectiveFunction.setParameters(allFreeParameters);
      cameraTworld = objectiveFunction.getPoseCameraTworld();
      intrinsics = objectiveFunction.getIntrinsics();
    }

    
    template <class FloatType, class Iter3D, class Iter2D>
    void
    estimateCameraParametersPinhole(
      CameraIntrinsicsPinhole<FloatType>& intrinsics,
      numeric::Transform3D<FloatType>& cameraTworld,
      unsigned int numPixelsX,
      unsigned int numPixelsY,
      Iter3D points3DBegin,
      Iter3D points3DEnd,
      Iter2D points2DBegin)
    {
      // find the best fit projection from world coordinates to image.
      //
      //        [[h0, h1, h2 ,  h3],
      //   H ==  [h4, h5, h6 ,  h7],
      //         [h8, h9, h10, h11]]
      //
      // xxx: Move estimateTransform3DTo2D into dlrNumeric? Make it
      // public somehow.
      numeric::Transform3DTo2D<FloatType> imagePworld =
        privateCode::estimateTransform3DTo2D<FloatType>(
          points3DBegin, points3DEnd, points2DBegin);
      FloatType h0 = imagePworld(0, 0);
      FloatType h1 = imagePworld(0, 1);
      FloatType h2 = imagePworld(0, 2);
      FloatType h3 = imagePworld(0, 3);
      FloatType h4 = imagePworld(1, 0);
      FloatType h5 = imagePworld(1, 1);
      FloatType h6 = imagePworld(1, 2);
      FloatType h7 = imagePworld(1, 3);
      FloatType h8 = imagePworld(2, 0);
      FloatType h9 = imagePworld(2, 1);
      FloatType h10 = imagePworld(2, 2);
      FloatType h11 = imagePworld(2, 3);
      
      // Assume H factors into:
      //
      //        [[f_u,   0, u_0, 0],   [[           ],
      //   H ~=  [  0, f_v, v_0, 0], *  [   R      t],
      //         [  0    0    1, 0]]    [           ],
      //                                [ 0, 0, 0, 1]]
      //
      // where "~=" denotes equality up to a scale factor, R is a
      // 3x3 rotation matrix, t is a 3x1 translation vector, and the
      // remaining variables are camera intrinsic parameters.
      //
      // Equivalently, we have:
      //
      //   [[a/f_u,     0,  -a*u_0/f_u],       [[           ],
      //    [    0, a/f_v,  -a*v_0/f_v], * H =  [   R      t],
      //    [    0,     0,           a]]        [           ]]
      //
      // where a is the unknown scale factor.  Inspired by Zhang, "A
      // Flexible New Technique for Camera Calibration," we rearrange
      // to get:
      //
      //         [[b0,  0, b1],
      //   J^T *  [ 0, b2, b3]  * J = R^T * R = I
      //          [b1, b3, b4]]
      //
      // where J is the first 3 columns of H, and b0...b4 are as follows:
      // 
      //   b0 = (a/f_u)^2
      //   b1 = -(a/f_u)^2 * u_0
      //   b2 = (a/f_v)^2
      //   b3 = -(a/f_v)^2 * v_0,
      //   b4 = (a/f_u)^2 * u_0^2 + (a/fv)^2 * v_0^2 + a^2
      //
      // Assuming non-singular J, this is easy to solve for b0 - b4:
      //
      //   [[b0,  0, b1],
      //    [ 0, b2, b3]  = J^(-T) *  J^(-1)
      //    [b1, b3, b4]]
      // 
      // We rearrange and solve linearly for b0 - b4 here.
      numeric::Array1D<FloatType> eyeVector(9);
      eyeVector = 0.0;
      eyeVector[0] = eyeVector[4] = eyeVector[8] = 1.0;

      numeric::Array2D<FloatType> coefficients(9, 5);
      coefficients(0, 0) = h0 * h0; coefficients(0, 1) = 2.0 * h0 * h8;
      coefficients(0, 2) = h4 * h4; coefficients(0, 3) = 2.0 * h4 * h8;
      coefficients(0, 4) = h8 * h8;
      coefficients(1, 0) = h0 * h1; coefficients(1, 1) = h8 * h1 + h0 * h9;
      coefficients(1, 2) = h4 * h5; coefficients(1, 3) = h8 * h5 + h4 * h9;
      coefficients(1, 4) = h8 * h9;
      coefficients(2, 0) = h0 * h2; coefficients(2, 1) = h8 * h2 + h0 * h10;
      coefficients(2, 2) = h4 * h6; coefficients(2, 3) = h8 * h6 + h4 * h10;
      coefficients(2, 4) = h8 * h10;
      coefficients(3, 0) = h0 * h1; coefficients(3, 1) = h0 * h9 + h8 * h1;
      coefficients(3, 2) = h4 * h5; coefficients(3, 3) = h4 * h9 + h8 * h5;
      coefficients(3, 4) = h8 * h9;
      coefficients(4, 0) = h1 * h1; coefficients(4, 1) = 2.0 * h1 * h9;
      coefficients(4, 2) = h5 * h5; coefficients(4, 3) = 2.0 * h5 * h9;
      coefficients(4, 4) = h9 * h9;
      coefficients(5, 0) = h1 * h2; coefficients(5, 1) = h9 * h2 + h1 * h10;
      coefficients(5, 2) = h5 * h6; coefficients(5, 3) = h9 * h6 + h5 * h10;
      coefficients(5, 4) = h9 * h10;
      coefficients(6, 0) = h0 * h2; coefficients(6, 1) = h0 * h10 + h8 * h2;
      coefficients(6, 2) = h4 * h6; coefficients(6, 3) = h4 * h10 + h8 * h6;
      coefficients(6, 4) = h8 * h10;
      coefficients(7, 0) = h1 * h2; coefficients(7, 1) = h1 * h10 + h9 * h2;
      coefficients(7, 2) = h5 * h6; coefficients(7, 3) = h5 * h10 + h9 * h6;
      coefficients(7, 4) = h9 * h10;
      coefficients(8, 0) = h2 * h2; coefficients(8, 1) = 2.0 * h2 * h10;
      coefficients(8, 2) = h6 * h6; coefficients(8, 3) = 2.0 * h6 * h10;
      coefficients(8, 4) = h10 * h10;

      // Solve for b0 - b4.
      numeric::Array1D<FloatType> bVector = 
        linearAlgebra::linearLeastSquares(coefficients, eyeVector);

      // Remembering the definitions of the elements of bVector, we
      // can (almost) solve algebraically for our camera parameters.
      // We can't solve for the sign of a, f_u, or f_v, since these
      // only show up in squared quanties.  These both have an impact
      // on the signs of the elements of the recovered transform, but
      // their effect on the total projection (through cameraTworld
      // and camera intrinsics) cancels out.  For this reason, we
      // assume two things: the focal length is positive; and the 3D
      // positions of the projected points are generally in front of
      // the camera.
      // 
      //   b0 = (a/f_u)^2
      //   b1 = -(a/f_u)^2 * u_0
      //   b2 = (a/f_v)^2
      //   b3 = -(a/f_v)^2 * v_0,
      //   b4 = (a/f_u)^2 * u_0^2 + (a/fv)^2 * v_0^2 + a^2
      if(bVector[0] < 0.0) {
        bVector *= -1.0;
      }
      if((bVector[0] < 1.0e-10) || (bVector[2] < 1.0e-10)) {
        BRICK_THROW(brick::common::ValueException,
                    "preEstimateCameraParameters()",
                    "Untenable initial estimate for intrinsics.");
      }
      FloatType u_0 = -bVector[1] / bVector[0];        
      FloatType v_0 = -bVector[3] / bVector[2];
      FloatType aSquared = bVector[4] + u_0 * bVector[1] + v_0 * bVector[3];
      if(aSquared < 1.0e-10) {
        BRICK_THROW(brick::common::ValueException,
                    "preEstimateCameraParameters()",
                    "Untenable scale in initial estimate for intrinsics.");
      }
      FloatType oneOverF_u = brick::common::squareRoot(bVector[0] / aSquared);
      // Note: assuming positive focal length.
      FloatType f_u = 1.0 / oneOverF_u;  
      FloatType oneOverF_v = brick::common::squareRoot(bVector[2] / aSquared);
      // Note: assuming positive focal length.
      FloatType f_v = 1.0 / oneOverF_v;  
      FloatType aa = brick::common::squareRoot(aSquared);
           
      // OK, I've copied one of the above equations here for easy
      // reference.  We just solved for the first matrix on the right
      // side of the equation, and we have H since the very beginning
      // of this function call.  This lets us solve for R and t.
      //
      //   [[a/f_u,     0,  -a*u_0/f_u],       [[           ],
      //    [    0, a/f_v,  -a*v_0/f_v], * H =  [   R      t],
      //    [    0,     0,           a]]        [           ]]
      //
      // This is a good time to notice the effect of having the signs
      // of a and/or f_u & f_v wrong.  Changing the sign of a simply
      // reverses the signs of all elements of R and t.  Changing the
      // sign of f_u reverses the sign of the first row of R and the
      // first element of t.  Changing the sign of f_v reverses the
      // sign of the second row of R and the second element of t.
      numeric::Array2D<FloatType> AInverse(3, 3);
      AInverse = 0.0;
      AInverse(0, 0) = aa / f_u;
      AInverse(0, 2) = -aa * u_0 / f_u;
      AInverse(1, 1) = aa / f_v;
      AInverse(1, 2) = -aa * v_0 / f_v;
      AInverse(2, 2) = aa;

      // We'll represent the 1st 3 columns of H separately.
      numeric::Array2D<FloatType> H0To3(3, 3);
      H0To3(0, 0) = h0; H0To3(0, 1) = h1; H0To3(0, 2) = h2;
      H0To3(1, 0) = h4; H0To3(1, 1) = h5; H0To3(1, 2) = h6;
      H0To3(2, 0) = h8; H0To3(2, 1) = h9; H0To3(2, 2) = h10;

      // Here's the final column of H.
      numeric::Array1D<FloatType> H3(3);
      H3[0] = h3; H3[1] = h7; H3[2] = h11;

      // Solve for the rigid transform.
      numeric::Array2D<FloatType> rotation =
        numeric::matrixMultiply<FloatType>(AInverse, H0To3);
      numeric::Array1D<FloatType> translation =
        numeric::matrixMultiply<FloatType>(AInverse, H3);

      // We know that the rotation matrix must be orthogonal (because
      // it's a rotation matrix).  It's pretty straightforward to
      // orthogonalalize using SVD.  This approach provably minimizes
      // the difference (usig Frobenius norm) between our original
      // estimate and the orthogonalized result.
      numeric::Array2D<FloatType> uMatrix;
      numeric::Array2D<FloatType> vTransposeMatrix;
      numeric::Array1D<FloatType> sigmaArray;
      linearAlgebra::singularValueDecomposition(
        rotation, uMatrix, sigmaArray, vTransposeMatrix);
      rotation = matrixMultiply<FloatType>(uMatrix, vTransposeMatrix);

      // Now we have our estimates all figured out, except for the
      // sign of scale parameter "a".  Go ahead and copy the result
      // into our output parameters.
      intrinsics = CameraIntrinsicsPinhole<FloatType>(
        numPixelsX, numPixelsY, 1.0, 1.0 / f_u, 1.0 / f_v, u_0, v_0);
      cameraTworld.setValue(
        rotation(0, 0), rotation(0, 1), rotation(0, 2), translation[0],
        rotation(1, 0), rotation(1, 1), rotation(1, 2), translation[1],
        rotation(2, 0), rotation(2, 1), rotation(2, 2), translation[2],
        0.0, 0.0, 0.0, 1.0);

      // Check to see whether, on average, our set of 3D points
      // projects from in front of the camera, or behind it.
      FloatType averageZValue =
        privateCode::computeAverageTransformedZValue(
          cameraTworld, points3DBegin, points3DEnd);
      
      if(averageZValue < 0.0) {
        // Based on our assumption of projecting from positive Z,
        // looks like we need to flip the sign of alpha.  We observed
        // earlier that this just changes the sign of the first three
        // rows of cameraTworld.
      cameraTworld.setValue(
        -rotation(0, 0), -rotation(0, 1), -rotation(0, 2), -translation[0],
        -rotation(1, 0), -rotation(1, 1), -rotation(1, 2), -translation[1],
        -rotation(2, 0), -rotation(2, 1), -rotation(2, 2), -translation[2],
        0.0, 0.0, 0.0, 1.0);
      }

      // Phew!
    }

  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CALIBRATIONTOOLS_IMPL_HH */
