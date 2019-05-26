#pragma once

#include "ProgramParams.h"

vector<complex<double>> DifferentialAlg(int N, Potential const & signal) { // N - number of nodes in the interval [T1,T2]
  double T1 = signal.T1();
  double T2 = signal.T2();
  double h = (T2 - T1) / N;
  complex<double> L = signal.L(T1);
  vector<complex<double>> B1; B1.reserve(N + 1);
  vector<complex<double>> B2; B2.reserve(N + 1);
  B1.push_back(0.);
  B2.push_back(conj(L));
  vector<complex<double>> q; q.reserve(N + 1);
  q.push_back(-2. * L);
  for (int j = 1; j <= N; j++) {
    vector<complex<double>> tmp_B1; tmp_B1.reserve(j + 1);
    vector<complex<double>> tmp_B2; tmp_B2.reserve(j + 1);
    tmp_B1.push_back(0.);
    tmp_B2.push_back(conj(L));
    for (int n = 1; n < j; n++) {
      tmp_B1.push_back(B1[n - 1] + h*q[j - 1] * B2[n - 1]);
      tmp_B2.push_back(B2[n] - h*conj(q[j - 1])*B1[n]);
    }
    B1 = tmp_B1;
    B2 = tmp_B2;
    B1[j] = B1[j - 1] + h*q[j - 1] * B2[j - 1];
    B2[j] = conj(signal.L(2. * (T1 + j*h)));
    for (int k = 0; k < j; k++) {
      B2[j] += conj(2.*h*B1[k] * signal.L(2.*(T1 + k*h)));
    }
    q.push_back(-2.*conj(B2[j]));
  }
  q.shrink_to_fit();
  return q;
}