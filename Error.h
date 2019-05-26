#pragma once
#include "ProgramParams.h"

void ErrorAnalyzer(const vector<complex<double>> & analytical,
                   const vector<complex<double>> & numerical, 
                   double T1, double h, const string & filepath) {
  assert(analytical.size() == numerical.size(), "sizes of vectors are not equal");
  ofstream file(filepath);
  for (int j = 0; j < numerical.size(); j++) {
    file << (T1 + j*h) << " " << abs((analytical[j] - numerical[j])) << endl;
  }
  file.close();
}

vector<double> ErrorAbs(const vector<complex<double>> & analytical,
                        const vector<complex<double>> & numerical  ) {
  size_t n = analytical.size();
  assert(n == numerical.size(), "sizes of vectors are not equal");
  vector<double> res(n);
  for (int j = 0; j < n; j++) {
    res[j] = abs(analytical[j] - numerical[j]);
  }
  return res;
}

double ErrorL1(const vector<complex<double>> & analytical,
  const vector<complex<double>> & numerical, const double h) {
  size_t n = analytical.size();
  assert(n == numerical.size(), "sizes of vectors are not equal");
  double res = 0.;
  for (int j = 0; j < n; j++) {
    res += abs(analytical[j] - numerical[j]);
  }
  res *= h;
  return res;
}

double ErrorL2(const vector<complex<double>> & analytical,
  const vector<complex<double>> & numerical, const double h) {
  size_t n = analytical.size();
  assert(n == numerical.size(), "sizes of vectors are not equal");
  double res = 0.;
  for (int j = 0; j < n; j++) {
    res += abs(analytical[j] - numerical[j])*abs(analytical[j] - numerical[j]);
  }
  res *= h;
  return res;
}