#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  //measurement matrix - laser (H_laser)
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0; 
  //measurement matrix - radar (Hj_) - does it really need to be initialized?
  /*
  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
         1, 1, 1, 1;
  */

  // Todo: Check if it needs to be set; for now I'm keeping it similar to the lecture notes with regards 
  // to which variables are initialized

  // state covariance matrix (see header)
  ekf_.P_ = MatrixXd(4, 4);
  
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;
  
  /*
  ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
  */
            
            
            

  // transition matrix
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;
  

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */

  // For initialization, either set the attributes directly or use the Init() function
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.

      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);

      // Todo: have to try different initialization methods, maybe for vx, vy?
      ekf_.x_ << rho * cos(phi), 
                 rho * sin(phi),
                 rho_dot*cos(phi),
                 rho_dot*sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      // Todo: Using the initialization from the lecture, could also try different ones for vx, vy
      ekf_.x_ << measurement_pack.raw_measurements_[0], 
              measurement_pack.raw_measurements_[1], 
              0, 
              0;

    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // updating state transition matrix F according to new elapsed time
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;  // elapsed time in seconds
  previous_timestamp_ = measurement_pack.timestamp_;  // update timestamp
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // updating process noise covariance matrix
  float noise_ax = 9;
  float noise_ay = 9;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << pow(dt, 4.0)/4.0 * noise_ax, 0, pow(dt, 3.0)/2.0 * noise_ax, 0,
            0, pow(dt, 4.0)/4.0 * noise_ay, 0, pow(dt, 3.0)/2.0 * noise_ay,
            pow(dt, 3.0)/2.0 * noise_ax, 0, pow(dt, 2.0) * noise_ax, 0,
            0, pow(dt, 3.0)/2.0 * noise_ay, 0, pow(dt, 2.0) * noise_ay;


  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    // set attributes of ekf_ for RADAR
    //Tools tools;  // initialize using constructor
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.H_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
    // set attributes of ekf_ for LASER
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
