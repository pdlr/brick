/**
***************************************************************************
* @file brick/computerVision/calibrationToolsRobust_impl.hh
*
* Header file defining inline and template functions declared in
* calibrationTools.hh.
*
* Copyright (C) 2010-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CALIBRATIONTOOLSROBUST_IMPL_HH
#define BRICK_COMPUTERVISION_CALIBRATIONTOOLSROBUST_IMPL_HH

// This file is included by calibrationToolsRobust.hh, and should not be
// directly included by user code, so no need to include
// calibrationToolsRobust.hh here.
//
// #include <brick/computerVision/calibrationToolsRobust.hh>


#include <brick/computerVision/calibrationTools.hh>
#include <brick/computerVision/cameraIntrinsicsPlumbBob.hh>
#include <brick/computerVision/ransacClassInterface.hh>

namespace brick {

  namespace computerVision {

    namespace privateCode {

      template <class FloatType>
      struct CTRSamplePair {
        brick::numeric::Vector3D<FloatType> point3D;
        brick::numeric::Vector2D<FloatType> point2D;

        CTRSamplePair()
          : point3D(), point2D() {}

        CTRSamplePair(brick::numeric::Vector3D<FloatType> a_point3D,
                      brick::numeric::Vector2D<FloatType> a_point2D)
          : point3D(a_point3D), point2D(a_point2D) {}
      };


      template <class Intrinsics>
      struct CTRResult {
        numeric::Transform3D<typename Intrinsics::FloatType> cameraTworld;
        Intrinsics intrinsics;
        CameraParameterEstimationStatistics<typename Intrinsics::FloatType>
          statistics;
      };


      // This abstract base class is the root of a tree that will be
      // used by RANSAC to do all of the problem-dependent work when
      // we call RANSAC from inside estimateCameraParameters*Robust().
      //
      // The first Ransac template argument says that we're going to
      // estimate the parameters based on a collection of
      // CTRSamplePair instances.  The second Ransac template argument
      // says that we're going to represent the result using a
      // CTRResult<Intrinsics> instance.
      template <class Intrinsics>
      class EstimateCameraParametersProblemParent
        : public RansacProblem< CTRSamplePair<typename Intrinsics::FloatType>,
                                CTRResult<Intrinsics> >
      {
      public:

        // Typedefs to keep syntax sane below.
        typedef CTRSamplePair<typename Intrinsics::FloatType> MySamplePair;
        typedef CTRResult<Intrinsics> MyResult;
        typedef RansacProblem<MySamplePair, MyResult> MyRansacProblem;
        typedef typename Intrinsics::FloatType FloatType;

        // SampleSequenceType is typedef'd inside RansacProblem.  It's
        // just a pair of iterators that define a sequence of
        // CTRSamplePair instances.
        typedef typename MyRansacProblem::SampleSequenceType
          MySampleSequenceType;


        // When constructing the EstimateCameraParametersProblemParent
        // instance, we have to pass in a sequence of CTRSamplePair
        // instances so that the RANSAC implementation can later ask
        // our EstimateCameraParametersProblemParent instance for
        // randomly selected sets of samples.  The first argument of
        // the parent class constructor, numSamplesRequired, specifies
        // how many samples (CTRSamplePair instances) are needed to
        // estimate the camera parameters.
        template <class IterType>
        EstimateCameraParametersProblemParent(
          Intrinsics const& intrinsics, FloatType maxResidual,
          unsigned int numSamplesRequired, IterType beginIter, IterType endIter)
          : MyRansacProblem(numSamplesRequired, beginIter, endIter),
            m_intrinsics(intrinsics),
            m_maxResidual(maxResidual),
            m_numPixelsX(intrinsics.getNumPixelsX()),
            m_numPixelsY(intrinsics.getNumPixelsY()) {}


        // This should be overridden by derived classes to estimate
        // camera parameters based on an input sequence of 3D-to-2D
        // point correspondences.
        virtual CTRResult<Intrinsics>
        estimateModel(MySampleSequenceType const& sampleSequence) = 0;


        // Given camera parameters, compute the residual for each
        // 3D-2D sample pair of the input sampleSequence.
        template <class IterType>
        void
        computeError(CTRResult<Intrinsics> const& model,
                     MySampleSequenceType const& sampleSequence,
                     IterType outputIter) {
          std::vector< brick::numeric::Vector3D<FloatType> > worldPoints;
          std::vector< brick::numeric::Vector2D<FloatType> > imagePoints;
          this->unpackSampleSequence(sampleSequence, worldPoints, imagePoints);

          for(unsigned int ii = 0; ii < worldPoints.size(); ++ii) {
            brick::numeric::Vector3D<FloatType> cameraPoint =
              model.cameraTworld * worldPoints[ii];
            brick::numeric::Vector2D<FloatType> projectedPoint =
              model.intrinsics.project(cameraPoint);
            *outputIter = brick::numeric::magnitude<FloatType>(
              projectedPoint - imagePoints[ii]);
            ++outputIter;
          }
        }


        // How closely (in pixels) must a world point project to its
        // corresponding image point in order to not be considered an
        // outlier.
        FloatType
        getNaiveErrorThreshold() {return m_maxResidual;}

      protected:

        void
        unpackSampleSequence(
          MySampleSequenceType const& sampleSequence,
          std::vector< brick::numeric::Vector3D<FloatType> >& worldPoints,
          std::vector< brick::numeric::Vector2D<FloatType> >& imagePoints) {

          // We have to copy sampleSequence so that we can increment the
          // "begin iterator" part of it.
          MySampleSequenceType mutableSequence = sampleSequence;

          unsigned int numberOfSamples =
            mutableSequence.second - mutableSequence.first;
          worldPoints.resize(numberOfSamples);
          imagePoints.resize(numberOfSamples);
          for(unsigned int ii = 0; ii < numberOfSamples; ++ii) {
            worldPoints[ii] = (mutableSequence.first)->point3D;
            imagePoints[ii] = (mutableSequence.first)->point2D;
            ++(mutableSequence.first);
          }
        }


        Intrinsics m_intrinsics;
        FloatType m_maxResidual;
        unsigned int m_numPixelsX;
        unsigned int m_numPixelsY;
      };


      // This class specializes EstimateCameraParametersProblemParent
      // for pinhole cameras.
      template <class Intrinsics>
      class EstimateCameraParametersProblem
        : public EstimateCameraParametersProblemParent<Intrinsics>
      {
      public:
        typedef typename Intrinsics::FloatType FloatType;
        typedef typename EstimateCameraParametersProblemParent<Intrinsics>
          ::MySampleSequenceType MyMySampleSequenceType;

        // See EstimateCameraParametersProblemParent for documentation.
        template <class IterType>
        EstimateCameraParametersProblem(
          Intrinsics const& intrinsics, FloatType maxResidual,
          unsigned int numSamplesRequired, IterType beginIter, IterType endIter)
          : EstimateCameraParametersProblemParent<Intrinsics>(
              intrinsics, maxResidual, numSamplesRequired,
              beginIter, endIter) {}


        // See EstimateCameraParametersProblemParent for documentation.
        CTRResult<Intrinsics>
        estimateModel(MyMySampleSequenceType const& sampleSequence) {
          std::vector< brick::numeric::Vector3D<FloatType> > worldPoints;
          std::vector< brick::numeric::Vector2D<FloatType> > imagePoints;
          this->unpackSampleSequence(sampleSequence, worldPoints, imagePoints);

          CTRResult<Intrinsics> result;
          result.intrinsics = this->m_intrinsics;
          estimateCameraParameters(
            result.intrinsics, result.cameraTworld,
            result.statistics,
            this->m_numPixelsX, this->m_numPixelsY,
            worldPoints.begin(), worldPoints.end(), imagePoints.begin());
          return result;
        }
      };


      // This class specializes EstimateCameraParametersProblemParent
      // for pinhole cameras.
      template <class FloatType>
      class EstimateCameraParametersPinholeProblem
        : public EstimateCameraParametersProblemParent<
            CameraIntrinsicsPinhole<FloatType> >
      {
      public:

        typedef CameraIntrinsicsPinhole<FloatType> MyIntrinsics;
        typedef EstimateCameraParametersProblemParent<MyIntrinsics>
          MyParentType;
        typedef typename MyParentType::MySampleSequenceType
          MySampleSequenceType;

        // See EstimateCameraParametersProblemParent for documentation.
        template <class IterType>
        EstimateCameraParametersPinholeProblem(
          MyIntrinsics const& intrinsics,
          FloatType maxResidual,
          IterType beginIter, IterType endIter)
          : EstimateCameraParametersProblemParent<MyIntrinsics>(
            intrinsics, maxResidual, 6, beginIter, endIter) {}


        // See EstimateCameraParametersProblemParent for documentation.
        CTRResult<MyIntrinsics>
        estimateModel(MySampleSequenceType const& sampleSequence) {
          std::vector< brick::numeric::Vector3D<FloatType> > worldPoints;
          std::vector< brick::numeric::Vector2D<FloatType> > imagePoints;
          this->unpackSampleSequence(sampleSequence, worldPoints, imagePoints);

          CTRResult<MyIntrinsics> result;
          estimateCameraParametersPinhole(
            result.intrinsics, result.cameraTworld,
            MyParentType::m_numPixelsX, MyParentType::m_numPixelsY,
            worldPoints.begin(), worldPoints.end(), imagePoints.begin());
          return result;
        }
      };


      template <class FloatType, class Iter3D, class Iter2D>
      std::vector< CTRSamplePair<FloatType> >
      buildSampleVector(Iter3D points3DBegin, Iter3D points3DEnd,
                        Iter2D points2DBegin)
      {
        unsigned int numberOfSamples = points3DEnd - points3DBegin;
        std::vector< CTRSamplePair<FloatType> > sampleVector(numberOfSamples);
        for(unsigned int ii = 0; ii < numberOfSamples; ++ii) {
          sampleVector[ii] = CTRSamplePair<FloatType>(
            *points3DBegin, *points2DBegin);
          ++points3DBegin;
          ++points2DBegin;
        }
        return sampleVector;
      }


      template <class Intrinsics, class Iter3D, class Iter2D>
      void
      estimateCameraParametersGeneralRobust(
        Intrinsics& intrinsics,
        brick::numeric::Transform3D<typename Intrinsics::FloatType>&
          cameraTworld,
        unsigned int numPixelsX,
        unsigned int numPixelsY,
        unsigned int numSamplesRequired,
        Iter3D points3DBegin,
        Iter3D points3DEnd,
        Iter2D points2DBegin,
        typename Intrinsics::FloatType maxResidual,
        typename Intrinsics::FloatType inlierProbability,
        typename Intrinsics::FloatType requiredConfidence,
        unsigned int minConsensusSetSize,
        unsigned int verbosity)
      {
        typedef typename Intrinsics::FloatType FloatType;
        std::vector< privateCode::CTRSamplePair<FloatType> > sampleVector =
          privateCode::buildSampleVector<FloatType>(
            points3DBegin, points3DEnd, points2DBegin);
        intrinsics.setNumPixelsX(numPixelsX);
        intrinsics.setNumPixelsY(numPixelsY);
        privateCode::EstimateCameraParametersProblem<Intrinsics> problem(
          intrinsics, maxResidual, numSamplesRequired,
          sampleVector.begin(), sampleVector.end());
        Ransac< privateCode::EstimateCameraParametersProblem<Intrinsics> >
          ransac(problem, minConsensusSetSize, requiredConfidence,
                 inlierProbability, verbosity);
        privateCode::CTRResult<Intrinsics> result = ransac.getResult();
        cameraTworld = result.cameraTworld;
        intrinsics = result.intrinsics;
      }

    } // namespace privateCode


    // The general case is not implemented.
    template <class Intrinsics, class Iter3D, class Iter2D>
    void
    estimateCameraParametersRobust(
      Intrinsics& /* intrinsics */,
      brick::numeric::Transform3D<typename Intrinsics::FloatType>&
        /* cameraTworld */,
      unsigned int /* numPixelsX */,
      unsigned int /* numPixelsY */,
      Iter3D /* points3DBegin */,
      Iter3D /* points3DEnd */,
      Iter2D /* points2DBegin */,
      typename Intrinsics::FloatType /* maxResidual */,
      typename Intrinsics::FloatType /* inlierProbability */,
      typename Intrinsics::FloatType /* requiredConfidence */,
      unsigned int /* minConsensusSetSize */,
      unsigned int /* verbosity */)
    {
      BRICK_THROW(common::NotImplementedException,
                "estimateCameraParametersRobust()",
                "This function template must be further specialized before you "
                "use it with this type, however specializing it is very easy "
                "to do.");
    }


    template <class Iter3D, class Iter2D>
    void
    estimateCameraParametersRobust(
      CameraIntrinsicsPlumbBob<double>& intrinsics,
      brick::numeric::Transform3D<double>& cameraTworld,
      unsigned int numPixelsX, unsigned int numPixelsY,
      Iter3D points3DBegin, Iter3D points3DEnd,
      Iter2D points2DBegin,
      double maxResidual,
      double inlierProbability,
      double requiredConfidence,
      unsigned int minConsensusSetSize = 0,
      unsigned int verbosity = 0)
    {
      privateCode::estimateCameraParametersGeneralRobust(
        intrinsics, cameraTworld, numPixelsX, numPixelsY, 8,
        points3DBegin, points3DEnd, points2DBegin, maxResidual,
        inlierProbability, requiredConfidence, minConsensusSetSize,
        verbosity);
    }


    template <class FloatType, class Iter3D, class Iter2D>
    void
    estimateCameraParametersPinholeRobust(
      CameraIntrinsicsPinhole<FloatType>& intrinsics,
      brick::numeric::Transform3D<FloatType>& cameraTworld,
      unsigned int numPixelsX,
      unsigned int numPixelsY,
      Iter3D points3DBegin,
      Iter3D points3DEnd,
      Iter2D points2DBegin,
      FloatType maxResidual,
      FloatType inlierProbability,
      FloatType requiredConfidence,
      unsigned int minConsensusSetSize,
      unsigned int verbosity)
    {
      std::vector< privateCode::CTRSamplePair<FloatType> > sampleVector =
        privateCode::buildSampleVector<FloatType>(
          points3DBegin, points3DEnd, points2DBegin);
      intrinsics.setNumPixelsX(numPixelsX);
      intrinsics.setNumPixelsY(numPixelsY);
      privateCode::EstimateCameraParametersPinholeProblem<FloatType> problem(
        intrinsics, maxResidual, sampleVector.begin(), sampleVector.end());
      Ransac< privateCode::EstimateCameraParametersPinholeProblem<FloatType> >
        ransac(
          problem, minConsensusSetSize, requiredConfidence, inlierProbability,
          verbosity);
      privateCode::CTRResult< CameraIntrinsicsPinhole<FloatType> > result =
        ransac.getResult();
      cameraTworld = result.cameraTworld;
      intrinsics = result.intrinsics;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CALIBRATIONTOOLSROBUST_IMPL_HH */
