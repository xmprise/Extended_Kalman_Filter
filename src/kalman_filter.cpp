#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

void KalmanFilter::Predict() {
  // x′ = Fx+u
  x_ = F_ * x_;
  // P′= FPF(transpose)+Q
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // y = z − Hx′
  VectorXd z_pred = (H_ * x_);
  VectorXd y = z - z_pred;

  Measurement(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  double px = x_[0];
  double py = x_[1];
  double vx = x_[2];
  double vy = x_[3];
  VectorXd z_pred = VectorXd(3);
  double sqrpxpy = pow(px, 2) + pow(py, 2);
  double sqrtpxpy = sqrt(sqrpxpy);

  if (px!=0 && sqrtpxpy > 0.001){
    z_pred << sqrtpxpy, atan2(py,px), (px*vx + py*vy)/sqrtpxpy;
  } else {
    z_pred << sqrtpxpy, 0, 0;
  }

  VectorXd y = z - z_pred;

  // normalize phi between -180 to 180
  y(1) = fmod(y(1)+180, 360);
  if(y(1) < 0)
    y(1)+=360;
  else
    y(1)-=180;

  Measurement(y);
}

void KalmanFilter::Measurement(const MatrixXd &y) {
  // S = HP′H(transpose)+R
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;

  // K = P′H(transpose)S(inverse)
  MatrixXd K = P_ * Ht * S.inverse();

  // x = x′ + Ky
  x_ = x_ + (K * y);

  // P=(I−KH)P′
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}