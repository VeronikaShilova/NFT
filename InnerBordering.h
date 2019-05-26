#pragma once

#include "ProgramParams.h"

vector<complex<double>> InnerBordering(int N, Potential const & signal) { // N - number of nodes in the interval [T1,T2]
  double tau1 = 2.*signal.T1();
  double tau2 = 2.*signal.T2();
  double h = (tau2 - tau1) / N;
  complex<double> L = signal.L(tau1);
  vector<complex<double>> y;  y.reserve(N + 1);
  y.push_back(1. / (1. + h*h*abs(L)*abs(L)*0.25));
  vector<complex<double>> z;  z.reserve(N + 1);
  z.push_back(-y[0] * h * L*0.5);
  vector<complex<double>> q; q.reserve(N + 1);
  q.push_back(-2. * L);

  for (int n = 1; n <= N; n++) {
    complex<double> b(0.);
    for (int j = 0; j < n; j++) {
      double tau = tau1 + (n - j)*h;
      b += h*signal.L(tau)*y[j];
    }
    q.push_back(-2.*b / h); // value of a signal at the moment (T1+n*h)
    double c = 1. / (1. + abs(b)*abs(b));
    complex<double> d = -b * c;
    y.push_back(0.);
    z.push_back(0.);
    vector<complex<double>> tmp_y = y;
    reverse(tmp_y.begin(), tmp_y.end());
    for_each(tmp_y.begin(), tmp_y.end(),
      [d](complex<double> &item) { item = conj(item); item *=  d; });
    vector<complex<double>> tmp_z = z;
    reverse(tmp_z.begin(), tmp_z.end());
    for_each(tmp_z.begin(), tmp_z.end(),
      [d](complex<double> &item) { item = conj(item); item *= -d; });
    for_each(y.begin(), y.end(),
      [c](complex<double> &item) { item *= c; });
    transform(tmp_z.begin(), tmp_z.end(), y.begin(), y.begin(), plus<complex<double>>());
    for_each(z.begin(), z.end(),
      [c](complex<double> &item) { item *= c; });
    transform(tmp_y.begin(), tmp_y.end(), z.begin(), z.begin(), plus<complex<double>>());
  }
  q.shrink_to_fit();
  return q;
}

