/**
***************************************************************************
* @file brick/computerVision/extendedKalmanFilter.hh
*
* Header file declaring an Extended Kalman Filter implementation.
*
* Copyright (C) 2009,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_EXTENDEDKALMANFILTER_HH
#define BRICK_COMPUTERVISION_EXTENDEDKALMANFILTER_HH

#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace computerVision {


    /**
     ** This class template implements the ExtendedKalman Filter[??].
     **
     ** This interface is not stable, and may not even stay as part of
     ** dlrComputerVision.  Caveat emptor.
     **
     **/
    template <class FloatType>
    class ExtendedKalmanFilter {
    public:

      /**
       * Default constructor.
       */
      explicit
      ExtendedKalmanFilter(FloatType startTime = 0.0);


      /**
       * Destructor.
       */
      virtual
      ~ExtendedKalmanFilter() {};


      /**
       * Use this member function tell the filter about a new
       * measurement, and to request that the state estimate be
       * updated to reflect this new measurement.  Under the hood, it
       * just calls doPredictionStep() and doMeasurementUpdate() in
       * sequence.
       *
       * @param measurementID This argument identifies to which
       * measurement model the measurment corresponds.
       *
       * @param timestamp This argument indicates the time at which
       * the measurement was acquired.
       *
       * @param measurement This argument specifies the value of the
       * measurement.
       *
       * @param controlInput This argument specifies the control input
       */
      virtual void
      addMeasurement(unsigned int measurementID,
                     FloatType timestamp,
                     brick::numeric::Array1D<FloatType> const& measurement,
                     brick::numeric::Array1D<FloatType> const& controlInput);


      /**
       * This member function provides a sanity check on the
       * user-provided Jacobian computation, which is a common source
       * of errors.  In its current version it only checks the
       * Jacobian with respect to the state vector.
       *
       * It works by computing a finite-difference approximation of
       * the Jacobian, and then finding the difference between this
       * approximation and the Jacobian provided by
       * this->getMeasurementJacobians().  The difference is returned
       * through a reference argument so that it can be evaluated by
       * the calling context.
       *
       * @param measurementID This argument identifies for which
       * measurement model the Jacobian should be checked.
       *
       * @param epsilonArray This argument specifies the step size for
       * the finite-difference approximation.  Each element of
       * epsilonArray corresponds to one dimension of the state
       * vector.  Finite differences will be calculated in dimension
       * "n" by evaluating the measurement model around the current
       * state at offsets of +/- epsilon[n].
       *
       * @param residualArray0 This argument returns the difference
       * between the appoximated Jacobian with respect to state, and
       * the corresponding Jacobian returned by
       * getMeasurementJacobians().
       */
      virtual void
      checkMeasurementJacobians(unsigned int measurementID,
                                brick::numeric::Array1D<FloatType> const& epsilonArray,
                                brick::numeric::Array2D<FloatType>& residualArray0);


      /**
       * This member function provides a sanity check on the
       * user-provided Jacobian computation, which is a common source
       * of errors.  In its current version it only checks the
       * Jacobian with respect to the state vector.
       *
       * It works by computing a finite-difference approximation of
       * the Jacobian, and then finding the difference between this
       * approximation and the Jacobian provided by
       * this->getProcessJacobians().  The difference is returned
       * through a reference argument so that it can be evaluated by
       * the calling context.
       *
       * @param epsilonArray This argument specifies the step size for
       * the finite-difference approximation.  Each element of
       * epsilonArray corresponds to one dimension of the state
       * vector.  Finite differences will be calculated in dimension
       * "n" by evaluating the measurement model around the current
       * state at offsets of +/- epsilon[n].
       *
       * @param controlInput This argument specifies what control
       * input should be used when calling applyProcessModel() to
       * generate the finite differences.
       *
       * @param residualArray0 This argument returns the difference
       * between the appoximated Jacobian with respect to state, and
       * the corresponding Jacobian returned by
       * getProcessJacobians().
       */
      virtual void
      checkProcessJacobians(brick::numeric::Array1D<FloatType> const& epsilonArray,
                            brick::numeric::Array1D<FloatType> const& controlInput,
                            brick::numeric::Array2D<FloatType>& residualArray);


      /**
       * Use this member function only if you are not using member
       * function addMeasurement().  It uses the current process model
       * to estimate the process state at the specified time.
       * Normally, this member function is called directly by
       * addMeasurement(), but it is exposed here so that the user has
       * more control.  For example, it may be useful to update the
       * filter state after calling doPredictionStep(), but before
       * calling doMeasurementUpdate().
       *
       * @param currentTime This argument indicates the time at which
       * the prediction should apply.  Traditionally (for a discrete
       * Kalman filter), the timestamp is incremented by 1 at each
       * step.
       *
       * @param controlInput This argument specifies the control input
       * in effect since the previous timestamp.
       */
      void
      doPredictionStep(FloatType currentTime,
                       brick::numeric::Array1D<FloatType> const& controlInput);


      /**
       * Use this member function only if you are not using member
       * function addMeasurement().  It causes the internal state to
       * be updated based on a new measurement.  Normally, this member
       * function is called directly by addMeasurement(), but it is
       * exposed here so that the user has more control.  For example,
       * it may be useful to update the filter state after calling
       * doPredictionStep(), but before calling doMeasurementUpdate().
       *
       * @param measurementID This argument identifies to which
       * measurement model the measurment corresponds.
       *
       * @param measurement This argument specifies the value of the
       * measurement.
       */
      void
      doMeasurementUpdate(unsigned int measurementID,
                          brick::numeric::Array1D<FloatType> const& measurement);


      /**
       * This member function may optionally be called to disable
       * updates of Kalman gain and estimation error covariance.  If
       * the process model matches the actual system well, the Kalman
       * gain and estimation error covariance will often converge to
       * nearly constant values.  In this case, disabling updates will
       * prevent unnecessary calculation.
       */
      virtual void
      freezeKalmanGain();


      /**
       * This member function returns the current state estimate for
       * the filter, as well as the estimated covariance of the state
       * estimate.
       *
       * @param timestamp This argument returns the time of the most
       * recent state estimate.
       *
       * @param state This argument returns the state estimate.
       *
       * @param covariance This argument returns the estimated
       * covariance.
       */
      virtual void
      getStateEstimate(FloatType& timestamp,
                       brick::numeric::Array1D<FloatType>& state,
                       brick::numeric::Array2D<FloatType>& covariance);


      /**
       * This member function sets the initial state estimate for the
       * filter, as well as the covariance of any Gaussian noise
       * reflected in the initial state estimate.
       *
       * @param timestamp This argument specifies the time of the
       * initial state estimate.
       *
       * @param state This argument specifies the initial state estimate.
       *
       * @param covariance This argument specifies the covariance
       * associated with the initial state estimate.
       */
      virtual void
      setStateEstimate(FloatType timestamp,
                       brick::numeric::Array1D<FloatType> const& state,
                       brick::numeric::Array2D<FloatType> const& covariance);


      /**
       * This member function may optionally be called to reverse the
       * effect of freezeKalmanGain().  That is, it re-enables updates
       * to the Kalman gain and estimation error covariance.  It does
       * not retroactively update the estimates, rather it simply
       * reenables the updates for subsequent calculations.
       */
      virtual void
      unfreezeKalmanGain();

    protected:

      /**
       * Subclasses should override this function to implement the
       * measurement model(s) relevant to the filter.  For example, a
       * (not extended) Kalman filter would make this function
       * calculate quantity(ies) having the form (H * currentState).
       *
       * @param measurementID This argument indicates which of the
       * (possibly many) measurement models to simulate.
       *
       * @param currentTime This argument indicates the time to which
       * state should be extrapolated.  In standard discrete KF and
       * EKF formulations, this argument is ignored, because the time
       * step is always 1.  We include it to accomodate time-dependent
       * process models.
       *
       * @param previousTime This argument indicates the time at
       * which the previous update occurred. In standard discrete KF and
       * EKF formulations, this argument is ignored, because the time
       * step is always 1.  We include it to accomodate time-dependent
       * process models.
       *
       * @param currentState This argument specifies the state
       * estimate at the current time.
       *
       * @return The return value is an estimate of relevant
       * measurements.
       */
      virtual
      brick::numeric::Array1D<FloatType>
      applyMeasurementModel(unsigned int measurementID,
                            FloatType currentTime,
                            FloatType previousTime,
                            brick::numeric::Array1D<FloatType> const& currentState) = 0;


      /**
       * Subclasses should override this function to implement the
       * process model to be tracked by the filter.  For example, a
       * (not extended) Kalman filter would make this function
       * calculate the quantity (A * previousState + B *
       * controlInput).
       *
       * @param currentTime This argument indicates the time to which
       * state should be extrapolated.  In standard discrete KF and
       * EKF formulations, this argument is ignored, because the time
       * step is always 1.  We include it to accomodate time-dependent
       * process models.
       *
       * @param previousTime This argument indicates the time at
       * which the previous update occurred. In standard discrete KF and
       * EKF formulations, this argument is ignored, because the time
       * step is always 1.  We include it to accomodate time-dependent
       * process models.
       *
       * @param previousState This argument specifies the state
       * estimate at the time of the most recent update.
       *
       * @param controlInput This argument specifies control input in
       * effect between the time of the last update and currentTime.
       *
       * @return The return value is an updated state estimate.
       */
      virtual
      brick::numeric::Array1D<FloatType>
      applyProcessModel(FloatType currentTime,
                        FloatType previousTime,
                        brick::numeric::Array1D<FloatType> const& previousState,
                        brick::numeric::Array1D<FloatType> const& controlInput) = 0;


      /**
       * This member function should be overridden by subclasses to
       * return the first derivatives of measurement model with respect
       * to the input state, and with respect to the measurement noise.
       *
       * @param measurementID This argument indicates for which of the
       * (possibly many) measurement models to compute Jacobians.
       *
       * @param currentTime The current time is provided to accomodate
       * time-dependent measurement models.
       *
       * @param previousTime The current time is provided to accomodate
       * time-dependent measurement models.
       *
       * @param state This argument specifies the current state
       * estimate, in case the measurement model is dependent on state
       * (watch for bad interactions if it is, since this moves us
       * even further from the original Discrete Kalman Filter).
       *
       * @param stateJacobian This argument is used to return a matrix
       * in which each row reflects the first derivatives of the
       * corresponding element of the measurement with respect to the
       * elements of the input state.
       *
       * @param noiseJacobian This argument is used to return a
       * matrix in which each row reflects the first derivatives of
       * the corresponding element of the output with respect to the
       * elements of the measurement noise.  This will often just be
       * an identity matrix.
       */
      virtual
      void
      getMeasurementJacobians(unsigned int measurementID,
                              FloatType currentTime,
                              FloatType previousTime,
                              brick::numeric::Array1D<FloatType> const& state,
                              brick::numeric::Array2D<FloatType>& stateJacobian,
                              brick::numeric::Array2D<FloatType>& noiseJacobian) = 0;


      /**
       * This member function should be overridden by subclasses to
       * return the first derivative of the process model with respect
       * to the input state, and with respect to the control input.
       *
       * @param currentTime The current time is provided to accomodate
       * time-varying measurement models.
       *
       * @param previousTime The current time is provided to accomodate
       * time-dependent measurement models.
       *
       * @param state This argument specifies the current state
       * estimate, in case the process model is dependent on state
       * (watch for bad interactions if it is, since this moves us
       * even further from the original Discrete Kalman Filter).
       *
       * @param stateJacobian This argument is used to return a square
       * matrix in which each row reflects the first derivatives of
       * the corresponding element of the output with respect to the
       * elements of the input state.
       *
       * @param noiseJacobian This argument is used to return a
       * not-necessarily-square matrix in which each row reflects the
       * first derivatives of the corresponding element of the output
       * with respect to the elements of the process noise.
       */
      virtual
      void
      getProcessJacobians(FloatType currentTime,
                          FloatType previousTime,
                          brick::numeric::Array1D<FloatType> const& state,
                          brick::numeric::Array2D<FloatType>& stateJacobian,
                          brick::numeric::Array2D<FloatType>& noiseJacobian) = 0;


      /**
       * This member function returns the (modeled) covariance of the
       * measurement noise.
       *
       * @param measurementID This argument indicates for which of the
       * (possibly many) measurement models to return the noise
       * covariance.
       *
       * @param currentTime The current time is provided to accomodate
       * time-dependent process models.
       *
       * @param previousTime The current time is provided to accomodate
       * time-dependent measurement models.
       *
       * @return The return value is the covariance matrix.
       */
      virtual
      brick::numeric::Array2D<FloatType>
      getMeasurementNoiseCovariance(unsigned int measurementID,
                                    FloatType currentTime,
                                    FloatType previousTime) = 0;


      /**
       * This member function returns the (modeled) covariance of the
       * process noise.
       *
       * @param currentTime The current time is provided to accomodate
       * time-dependent process models.
       *
       * @param previousTime The current time is provided to accomodate
       * time-dependent measurement models.
       *
       * @return The return value is the covariance matrix.
       */
      virtual
      brick::numeric::Array2D<FloatType>
      getProcessNoiseCovariance(FloatType currentTime, FloatType previousTime) = 0;


      /* =================== Member variables =================== */

      brick::numeric::Array2D<FloatType> m_covariance;
      brick::numeric::Array1D<FloatType> m_state;
      FloatType m_previousTimestamp;
      FloatType m_timestamp;
    };


  } // namespace computerVision

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/extendedKalmanFilter_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_EXTENDEDKALMANFILTER_HH */
