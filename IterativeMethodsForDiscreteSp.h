#pragma once

#include "ProgramParams.h"

complex<double> Secant(Potential const & signal, complex<double> xi, int n) {
  complex<double> z, z1 = xi + 1e-05;
  vector<complex<double>> a1(2);
  vector<complex<double>> a2(2);
  while (abs(z - xi) > 1e-10) {
    a1 = BOsborne(signal, xi, n);
    a2 = BOsborne(signal, z1, n);
    z = z1 - a2[0] * ((z1 - xi) / (a2[0] - a1[0]));
    xi = z1;
    z1 = z;
  }
  return z;
}
