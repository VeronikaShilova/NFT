#pragma once

#include "ProgramParams.h"

class RectangularContour {
protected:
  const double Xmin_;
  const double Xmax_;
  const double Ymin_;
  const double Ymax_;
  int n_;
  double h_;
  int iSide1, iSide2, iSide3;
public:
  RectangularContour(const double Xmin, const double Xmax,
    const double Ymin, const double Ymax, int n)
    : Xmin_(Xmin), Xmax_(Xmax)
    , Ymin_(Ymin), Ymax_(Ymax) {
    setPointNumber(n);
  }
  void setPointNumber(int n) {
    n_ = n;
    h_ = 2. * (Xmax_ - Xmin_ + Ymax_ - Ymin_) / n_;
    double ih_ = 1. / h_;
    iSide1 = static_cast<int>(ih_ * (Xmax_ - Xmin_));
    iSide2 = static_cast<int>(ih_ * (Ymax_ - Ymin_));
    iSide3 = iSide1;
  }
  int getPointNumber() {
    return n_;
  }
  complex<double> z(int m) const {
    assert(m < n_);
    if ((m -= iSide1) <= 0) return complex<double>(Xmin_ + (m + iSide1 + 0.5) * h_, Ymin_); else
      if ((m -= iSide2) <= 0) return complex<double>(Xmax_, Ymin_ + (m + iSide2 + 0.5) * h_); else
        if ((m -= iSide3) <= 0) return complex<double>(Xmax_ - (m + iSide3 + 0.5) * h_, Ymax_); else
          return complex<double>(Xmin_, Ymax_ - (m + 0.5) * h_);
  }
};

/*class RectangularContour {
protected:
  const double Xmin_;
  const double Xmax_;
  const double Ymin_;
  const double Ymax_;
public:
  RectangularContour(const double Xmin, const double Xmax, 
                     const double Ymin, const double Ymax) : Xmin_(Xmin), Xmax_(Xmax), 
                                                             Ymin_(Ymin), Ymax_(Ymax) {}
  vector<complex<double>> getPoints(const int n) const {
    vector<complex<double>> z(n);
    const double h = ((Xmax_ - Xmin_) * 2 + (Ymax_ - Ymin_) * 2) / n;
    int j = 0, k = 1, m = 0;
    while (Xmin_ + j*h <= Xmax_) {
      z[m] = Xmin_ + j*h + i*Ymin_;
      m++;
      j++;
    }
    while (Ymin_ + k*h < Ymax_) {
      z[m] = Xmax_ + i*(Ymin_ + k*h);
      m++;
      k++;
    }
    j = 0;
    while (Xmax_ - j*h >= Xmin_) {
      z[m] = Xmax_ - j*h + i*Ymax_;
      m++;
      j++;
    }
    k = 1;
    while (Ymax_ - k*h > Ymin_) {
      z[m] = Xmin_ + i*(Ymax_ - k*h);
      m++;
      k++;
    }
    return z;
  }
};*/
