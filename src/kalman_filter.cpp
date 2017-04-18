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

    double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
    double phi = atan2(x_(1), x_(0));
    double rhodot;

    if (fabs(rho) < 0.0001) {
        rhodot = 0;
    } else {
        rhodot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
    }

    VectorXd z_pred = VectorXd(3);
    z_pred << rho, phi, rhodot;
    VectorXd y = z - z_pred;

    Measurement(y);
}

void KalmanFilter::Measurement(const MatrixXd &y) {
    // S = HP′H(transpose)+R
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;

    // K = P′H(transpose)S(inverse)
    MatrixXd K = PHt * Si;

    // x = x′ + Ky
    x_ = x_ + (K * y);

    // P=(I−KH)P′
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}