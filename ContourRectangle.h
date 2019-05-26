#pragma once

#include "ProgramParams.h"

class RectangularContour {
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
};
