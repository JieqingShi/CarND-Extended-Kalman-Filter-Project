#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size()==0){
      std::cout<<"CalculateRMSE() SizeZeroError!"<<std::endl;
      return rmse;
  }
  
  if(estimations.size()!=ground_truth.size()){
      std::cout<<"CalculateRMSE() NotEqualSizeError!"<<std::endl;
      return rmse;
  }

  for (int i=0; i < estimations.size(); ++i) {
    VectorXd error = estimations[i] - ground_truth[i];
    error = error.array() * error.array();  //this is how to square arrays!
    rmse += error;
  }

  rmse = rmse/estimations.size();

  rmse = rmse.array().sqrt();

  return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float denom = pow(px, 2.0)+pow(py, 2.0);
  float vxpyMinusvypx = vx * py - vy * px;
  if(fabs(denom)<0.0001){
      std::cout << "CalculateJacobian () - DivByZeroError!" << std::endl;
      return Hj;
  }
  Hj << px/sqrt(denom), py/sqrt(denom), 0, 0,
        -py/denom, px/denom, 0, 0,
        py*vxpyMinusvypx/pow(denom, 3.0/2.0), -px*vxpyMinusvypx/pow(denom, 3.0/2.0), px/sqrt(denom), py/sqrt(denom);
  return Hj;
}
