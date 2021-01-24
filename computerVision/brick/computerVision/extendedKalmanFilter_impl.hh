/**
***************************************************************************
* @file brick/computerVision/extendedKalmanFilter.cc
*
* Header file defining inline and template functions declared in
* extendedKalmanFilter.hh.
*
* Copyright (C) 2009,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_EXTENDEDKALMANFILTER_IMPL_HH
#define BRICK_COMPUTERVISION_EXTENDEDKALMANFILTER_IMPL_HH

// This file is included by extendedKalmanFilter.hh, and should not be
// directly included by user code, so no need to include
// extendedKalmanFilter.hh here.
//
// #include <brick/computerVision/extendedKalmanFilter.hh>

#include <brick/computerVision/extendedKalmanFilter.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace computerVision {

    // Default constructor.
    template <class FloatType>
    ExtendedKalmanFilter<FloatType>::
    ExtendedKalmanFilter(FloatType startTime)
      : m_covariance(),
        m_state(),
        m_previousTimestamp(startTime),
        m_timestamp(startTime)
    {
      // Empty.
    }


    // Use this member function tell the filter about a new
    // measurement, and to request that the state estimate be
    // updated to reflect this new measurement.
    template <class FloatType>
    void
    ExtendedKalmanFilter<FloatType>::
    addMeasurement(unsigned int measurementID,
                   FloatType timestamp,
                   brick::numeric::Array1D<FloatType> const& measurement,
                   brick::numeric::Array1D<FloatType> const& controlInput)
    {
      this->doPredictionStep(timestamp, controlInput);
      this->doMeasurementUpdate(measurementID, measurement);
    }


    template <class FloatType>
    void
    ExtendedKalmanFilter<FloatType>::
    checkMeasurementJacobians(unsigned int measurementID,
                              brick::numeric::Array1D<FloatType> const& epsilonArray,
                              brick::numeric::Array2D<FloatType>& residualArray0)
    {
      // Check arguments.
      if(epsilonArray.size() != m_state.size()) {
        BRICK_THROW(common::ValueException,
                  "ExtendedKalmanFilter::checkProcessJacobians()",
                  "Argument epsilon has incorrect size.");
      }

      // Get the symbolically computed Jacobians and check their sizes.
      brick::numeric::Array2D<FloatType> measurementJacobian0;
      brick::numeric::Array2D<FloatType> measurementJacobian1;
      this->getMeasurementJacobians(
        measurementID, m_timestamp, m_previousTimestamp, m_state,
        measurementJacobian0, measurementJacobian1);


      brick::numeric::Array1D<FloatType> measurement =
        this->applyMeasurementModel(
          measurementID, m_timestamp, m_previousTimestamp, m_state);
      brick::numeric::Array2D<FloatType> measurementNoiseCovariance =
        this->getMeasurementNoiseCovariance(
          measurementID, m_timestamp, m_previousTimestamp);
      if(measurementJacobian0.rows() != measurement.size()
         || measurementJacobian0.columns() != m_state.size()
         || measurementJacobian1.rows() != measurement.size()
         || (measurementJacobian1.columns()
             != measurementNoiseCovariance.rows())) {
        BRICK_THROW(common::LogicException,
                  "ExtendedKalmanFilter::checkMeasurementJacobians()",
                  "At least one of the returned Jacobians has incorrect size.");
      }

      // Generate finite-difference approximation to measurementJacobian0, the
      // Jacobian wrt state.
      brick::numeric::Array2D<FloatType> approximateJacobian0(
        measurementJacobian0.rows(), measurementJacobian0.columns());
      for(unsigned int columnIndex = 0; columnIndex < m_state.size();
          ++columnIndex) {
        brick::numeric::Array1D<FloatType> state = m_state.copy();

        state[columnIndex] -= epsilonArray[columnIndex];
        brick::numeric::Array1D<FloatType> resultMeasurement0 =
          this->applyMeasurementModel(
            measurementID, m_timestamp, m_previousTimestamp, state);

        state[columnIndex] += 2.0 * epsilonArray[columnIndex];
        brick::numeric::Array1D<FloatType> resultMeasurement1 =
          this->applyMeasurementModel(
            measurementID, m_timestamp, m_previousTimestamp, state);

        brick::numeric::Array1D<FloatType> approximateGradient =
          ((resultMeasurement1 - resultMeasurement0)
           / (2.0 * epsilonArray[columnIndex]));
        for(unsigned int rowIndex = 0; rowIndex < resultMeasurement0.size();
            ++rowIndex) {
          approximateJacobian0(rowIndex, columnIndex) =
            approximateGradient[rowIndex];
        }
      }

      // Unfortunately, we don't have a good way to check the Jacobian
      // wrt noise right now.  On the bright side, it's often the
      // identity matrix, which is hard to get wrong.  In cases where
      // it's not the identity matrix, we'll trust the user to write
      // her own test.

      residualArray0 = approximateJacobian0 - measurementJacobian0;
    }


    template <class FloatType>
    void
    ExtendedKalmanFilter<FloatType>::
    checkProcessJacobians(brick::numeric::Array1D<FloatType> const& epsilonArray,
                          brick::numeric::Array1D<FloatType> const& controlInput,
                          brick::numeric::Array2D<FloatType>& residualArray0)
    {
      // Check arguments.
      if(epsilonArray.size() != m_state.size()) {
        BRICK_THROW(common::ValueException,
                  "ExtendedKalmanFilter::checkProcessJacobians()",
                  "Argument epsilon has incorrect size.");
      }

      // Get the symbolically computed Jacobians and check their sizes.
      brick::numeric::Array2D<FloatType> processJacobian0;
      brick::numeric::Array2D<FloatType> processJacobian1;
      this->getProcessJacobians(
        m_timestamp, m_previousTimestamp, m_state,
        processJacobian0, processJacobian1);

      brick::numeric::Array2D<FloatType> processNoiseCovariance =
        this->getProcessNoiseCovariance(m_timestamp, m_previousTimestamp);
      if(processJacobian0.rows() != m_state.size()
         || processJacobian0.columns() != m_state.size()
         || processJacobian1.rows() != m_state.size()
         || processJacobian1.columns() != processNoiseCovariance.rows()) {
        BRICK_THROW(common::LogicException,
                  "ExtendedKalmanFilter::checkProcessJacobians()",
                  "At least one of the returned Jacobians has incorrect size.");
      }

      // Generate finite-difference approximation to processJacobian0, the
      // Jacobian wrt state.
      brick::numeric::Array2D<FloatType> approximateJacobian0(
        processJacobian0.rows(), processJacobian0.columns());
      for(unsigned int columnIndex = 0; columnIndex < m_state.size();
          ++columnIndex) {
        brick::numeric::Array1D<FloatType> state = m_state.copy();

        state[columnIndex] -= epsilonArray[columnIndex];
        brick::numeric::Array1D<FloatType> resultState0 =
          this->applyProcessModel(m_timestamp, m_previousTimestamp, state,
                                  controlInput);

        state[columnIndex] += 2.0 * epsilonArray[columnIndex];
        brick::numeric::Array1D<FloatType> resultState1 =
          this->applyProcessModel(m_timestamp, m_previousTimestamp, state,
                                  controlInput);

        brick::numeric::Array1D<FloatType> approximateGradient =
          (resultState1 - resultState0) / (2.0 * epsilonArray[columnIndex]);
        for(unsigned int rowIndex = 0; rowIndex < m_state.size(); ++rowIndex) {
          approximateJacobian0(rowIndex, columnIndex) =
            approximateGradient[rowIndex];
        }
      }

      // Unfortunately, we don't have a good way to check the Jacobian
      // wrt noise right now.  On the bright side, it's often the
      // identity matrix, which is hard to get wrong.  In cases where
      // it's not the identity matrix, we'll trust the user to write
      // her own test.

      residualArray0 = approximateJacobian0 - processJacobian0;
    }


    // This member function may optionally be called to disable
    // updates of Kalman gain and estimation error covariance.
    template <class FloatType>
    void
    ExtendedKalmanFilter<FloatType>::
    freezeKalmanGain()
    {
      BRICK_THROW(common::NotImplementedException,
                "ExtendedKalmanFilter::freezeKalmanGain()",
                "Sorry. Should be easy to implement, though.");
    }


    // This member function returns the current state estimate for
    // the filter, as well as the estimated covariance of the state
    // estimate.
    template <class FloatType>
    void
    ExtendedKalmanFilter<FloatType>::
    getStateEstimate(FloatType& timestamp,
                     brick::numeric::Array1D<FloatType>& state,
                     brick::numeric::Array2D<FloatType>& covariance)
    {
      timestamp = m_timestamp;
      state = m_state.copy();
      covariance = m_covariance.copy();
    }


    // This member function sets the initial state estimate for the
    // filter, as well as the covariance of any Gaussian noise
    // reflected in the initial state estimate.
    template <class FloatType>
    void
    ExtendedKalmanFilter<FloatType>::
    setStateEstimate(FloatType timestamp,
                     brick::numeric::Array1D<FloatType> const& state,
                     brick::numeric::Array2D<FloatType> const& covariance)
    {
      m_timestamp = timestamp;
      m_state = state.copy();
      m_covariance = covariance.copy();
    }


    // This member function may optionally be called to reverse the
    // effect of freezeKalmanGain().
    template <class FloatType>
    void
    ExtendedKalmanFilter<FloatType>::
    unfreezeKalmanGain()
    {
      BRICK_THROW(common::NotImplementedException,
                "ExtendedKalmanFilter::unfreezeKalmanGain()",
                "Sorry. Should be easy to implement, though.");
    }


    template <class FloatType>
    void
    ExtendedKalmanFilter<FloatType>::
    doPredictionStep(FloatType currentTime,
                     brick::numeric::Array1D<FloatType> const& controlInput)
    {
      // Prediction step is as follows:
      //
      //   xpr_k = f(xpo_(k-1), u_(k-1), 0)
      //
      // where xpr_k is the a priori (before adding in measurements
      // taken at time k) estimate of x at time k.  It will be updated
      // by doMeasurementUpdate() to generate the a posteriori
      // estimate at time k, which reflects any new sensor data.
      // xpo_(k-1) is the a posteriori estimate of x at time (k - 1),
      // u_(k-1) is the control input at time (k - 1), f() is the user
      // supplied processmodel (which can be nonlinear), and the final
      // argument of zero reflects that the prediction step assumes no
      // noise.
      //
      // Prediction step updates state covariance estimate as follows:
      //
      //   Ppr_k = A_k * Ppo_(k-1) * (A_k)^T + W_k * Q_(k-1) * (W_k)^T
      //
      // Where Ppr_k is the a priori estimate of the covariance of our
      // state estimate at time k, Ppo_(k-1) is the a posteriori
      // estimate of the covariance at time (k-1), Q_(k-1) is the
      // covariance of our process noise at time (k - 1), A_k and W_k
      // are the process Jacobians at time k, and the right
      // superscript ^T indicates matrix transpose.
      m_previousTimestamp = m_timestamp;
      m_timestamp = currentTime;
      m_state = this->applyProcessModel(
        m_timestamp, m_previousTimestamp, m_state, controlInput);

      // Note(xxx): Inefficient to transpose every time.
      brick::numeric::Array2D<FloatType> processJacobian0;
      brick::numeric::Array2D<FloatType> processJacobian1;
      this->getProcessJacobians(
        m_timestamp, m_previousTimestamp, m_state,
        processJacobian0, processJacobian1);
      m_covariance =
        brick::numeric::matrixMultiply<FloatType>(
          brick::numeric::matrixMultiply<FloatType>(
            processJacobian0, m_covariance),
          processJacobian0.transpose());
      m_covariance +=
        brick::numeric::matrixMultiply<FloatType>(
          brick::numeric::matrixMultiply<FloatType>(
            processJacobian1,
            this->getProcessNoiseCovariance(m_timestamp, m_previousTimestamp)),
          processJacobian1.transpose());
    }


    template <class FloatType>
    void
    ExtendedKalmanFilter<FloatType>::
    doMeasurementUpdate(unsigned int measurementID,
                        brick::numeric::Array1D<FloatType> const& measurement)
    {
      // Measurement step is as follows (assuming prediction step has
      // already brought us to the time of the sensor observation):
      //
      //   K_k = Ppr_k * (H_k)^T * (H_k * Ppr_k * (H_k)^T
      //                            + V_k * R * (V_k)^T)^(-1)
      //
      //   xpo_k = xpr_k + K_k * (z_k - H * xpr_k)
      //
      //   Ppo_k = (1 - K_k * H)Ppr_k
      //
      // Where K_k is the "Kalman gain" at time k, Ppr_k is the a
      // priori estimate of the covariance of our state estimate at
      // time k, H^T is the transpose of measurement matrix H, R is
      // the covariance of the measurement noise, xpr_k is the a
      // priori (before adding in measurements taken at time k)
      // estimate of x at time k, xpo_k is the a posteriori (after
      // accounting for measurements) estimate of x at time k, z_k is
      // the measurement at time k, and Ppo_k is the a posteriori
      // estimate of the covariance of our state estimate at time k.

      brick::numeric::Array2D<FloatType> measurementJacobian0;
      brick::numeric::Array2D<FloatType> measurementJacobian1;
      this->getMeasurementJacobians(
        measurementID, m_timestamp, m_previousTimestamp, m_state,
        measurementJacobian0, measurementJacobian1);

      // Convenient matrices for constructing Kalman gain, etc.
      brick::numeric::Array2D<FloatType> const& HMatrix = measurementJacobian0;
      brick::numeric::Array2D<FloatType> HTranspose = HMatrix.transpose();
      brick::numeric::Array2D<FloatType> const& VMatrix = measurementJacobian1;
      brick::numeric::Array2D<FloatType> VTranspose = VMatrix.transpose();

      // "Denominator" should be (H_k * Ppr_k * (H_k)^T) + (V_k * R_k * (V_k)^T)
      brick::numeric::Array2D<FloatType> denominator =
        brick::numeric::matrixMultiply<FloatType>(
          brick::numeric::matrixMultiply<FloatType>(HMatrix, m_covariance), HTranspose);
      denominator +=
        (brick::numeric::matrixMultiply<FloatType>(
          brick::numeric::matrixMultiply<FloatType>(
            VMatrix,
            this->getMeasurementNoiseCovariance(
              measurementID, m_timestamp, m_previousTimestamp)),
          VTranspose));

      brick::numeric::Array1D<FloatType> innovation =
        measurement - this->applyMeasurementModel(
          measurementID, m_timestamp, m_previousTimestamp, m_state);

      brick::numeric::Array2D<FloatType> kalmanGain =
        brick::numeric::matrixMultiply<FloatType>(
          brick::numeric::matrixMultiply<FloatType>(m_covariance, HTranspose),
                       linearAlgebra::inverse(denominator));

      m_state += brick::numeric::matrixMultiply<FloatType>(kalmanGain, innovation);

      // Compute P_k = (I - K * H) * P_(k-1)
      brick::numeric::Array2D<FloatType> kTimesHMinusI =
        brick::numeric::matrixMultiply<FloatType>(kalmanGain, HMatrix);
      if(kTimesHMinusI.rows() != kTimesHMinusI.columns()) {
        BRICK_THROW(common::LogicException,
                  "ExtendedKalmanFilter::doMeasurementUpdate()",
                  "Programming error. K * H must be square.");
      }
      for(unsigned int ii = 0; ii < kTimesHMinusI.rows(); ++ii) {
        kTimesHMinusI(ii, ii) -= 1.0;
      }
      m_covariance = brick::numeric::matrixMultiply<FloatType>(kTimesHMinusI, m_covariance);
      m_covariance *= -1.0;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_EXTENDEDKALMANFILTER_IMPL_HH */
