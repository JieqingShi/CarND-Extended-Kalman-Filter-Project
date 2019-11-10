#include <math.h>
#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

// have to make sure that all attributes x_, P_, F_, ... are set
void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;  // assuming u is always VectorXd(0,0);
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y = z - H_ * x_;
  UpdateMatrices(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float rho = sqrt(pow(px, 2.0) + pow(py, 2.0));
  float phi = atan2(py, px);
  
  float rho_dot = (px*vx + py*vy) / rho;
  VectorXd hx(3);
  hx << rho, phi, rho_dot;

  VectorXd y = z - hx;

  // making sure that the resulting angle phi in y is between -M_PI and M_PI
  while (y(1) > M_PI) {
      y(1) -= 2*M_PI;
  }
  while (y(1) < -M_PI) {
      y(1) += 2*M_PI;
  }
  
  UpdateMatrices(y);
}

void KalmanFilter::UpdateMatrices(const VectorXd &y) {
  /**
   * Update of matrices in the Update() and UpdateEKF() functions;
   * Outsourced the code into this function to reduce amount of Copy and Paste
   */

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  x_ = x_ + (K*y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H_) * P_; 
}
