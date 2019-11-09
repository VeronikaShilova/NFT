#pragma once

#include "ProgramParams.h"

//TODO: SHIFT TO EIGEN
vector<complex<double>> BOsborne(Potential const & signal, const complex<double> xi, const int n) {
  vector<complex<double>> ab(2);
  const double T1 = signal.T1();
  const double T2 = signal.T2();
  const double dt = (T2 - T1) / n;
  complex<double> f1(exp(-i * xi * T1)), f2(0.);

  for (int j = 1; j <= n; j++) {
    double t = T1 + (j - 0.5) * dt;
    complex<double> q = signal.q(t);
    double abs_q = abs(q);
    complex<double> tmp1(f1), tmp2(f2); 
    complex<double> k = sqrt(-abs_q*abs_q - xi*xi);
    f1 = (cosh(k * dt) - (i * xi / k) * sinh(k * dt)) * tmp1 + tmp2 * sinh(k * dt) * q / k;
    f2 = -tmp1 * conj(q) * sinh(k * dt) / k + (cosh(k * dt) + (i * xi / k) * sinh(k * dt)) * tmp2;
  }

  ab[0] = f1 * exp(i * xi * T2); // coefficient a
  ab[1] = f2 * exp(-i * xi * T2); // coefficient b
  return ab;
}

complex<double> BOsbornePrime(Potential const & signal, const complex<double> xi, const int n) {
  const double T1 = signal.T1();
  const double T2 = signal.T2();
  const double dt = (T2 - T1) / n;
  complex<double> f1(exp(-i * xi * T1)), f2(0.), fPrime1(-i * T1 * exp(-i * xi * T1)), fPrime2(0.);

  for (int j = 1; j <= n; j++) {
    double t = T1 + (j - 0.5) * dt;
    complex<double> q = signal.q(t);
    double abs_q = abs(q);
    complex<double> tmp1(f1), tmp2(f2), tmpPrime1(fPrime1), tmpPrime2(fPrime2);
    complex<double> k = sqrt(-abs_q * abs_q - xi * xi);
    complex<double> ik = i / k;
    complex<double> tPrime11(i * xi * xi * dt * cosh(k * dt) / (k * k) -
        (xi * dt / k + ik + i * xi * xi / (k * k * k)) * sinh(k * dt)),
                    tPrime12(xi * sinh(k * dt) * q / (k * k * k) - xi * dt * cosh(k * dt) * q / (k * k)),
                    tPrime21(-xi * conj(q) * sinh(k * dt) / (k * k * k) +
        xi * dt * conj(q) * cosh(k * dt) / (k * k)),
                    tPrime22(-i * xi * xi * dt * cosh(k * dt) / (k * k) +
        (-xi * dt / k + ik + i * xi * xi / (k * k * k)) * sinh(k * dt));
    complex<double> t11(cosh(k * dt) - (i * xi / k) * sinh(k * dt)),
                    t12(sinh(k * dt) * q / k), 
                    t21(-conj(q) * sinh(k * dt) / k),
                    t22(cosh(k * dt) + (i * xi / k) * sinh(k * dt));
    fPrime1 = tPrime11 * tmp1 + tPrime12 * tmp2 + t11 * tmpPrime1 + t12 * tmpPrime2;
    fPrime2 = tPrime21 * tmp1 + tPrime22 * tmp2 + t21 * tmpPrime1 + t22 * tmpPrime2;
    f1 = t11 * tmp1 + t12 * tmp2;
    f2 = t21 * tmp1 + t22 * tmp2;
  }

  complex<double> aPrime;
  aPrime = (fPrime1 + i * T2 * f1) * exp(i * xi * T2);
  return aPrime;
}
