#pragma once

#include "ProgramParams.h"

const double epsilon = 1e-1;
const double accuracy = 1e-10;

// intial values for iteration process are uniformly distributed on a circle with radius (rootMax-epsilon)
vector<complex<double>> Roots(int n, double rootMax) { 
  vector<complex<double>> roots(n);
  double L = 2.*pi*(rootMax - epsilon);
  double h = L / n;
  double phi;
  for (int j = 0; j < n; j++) {
    phi = (h*j) / (rootMax - epsilon);
    roots[j] = (rootMax - epsilon) * exp(i * phi);
  }
  return roots;
}

complex<double> Poly(int n, vector<complex<double>> sigma, complex<double> z) {
  complex<double> poly = pow(z, n);
  for (int j = 0; j < n; j++) {
    poly += pow(z, n - j - 1)*sigma[j];
  }
  return poly;
}

vector<complex<double>> findRoots(int n, vector<complex<double>> sigma) {
  double rootMax = abs(*max_element(begin(sigma), end(sigma), 
    [](complex<double> item1, complex<double> item2) {return abs(item1) < abs(item2); })) + 1.;
  vector<complex<double>> roots = Roots(n, rootMax);
  complex<double> denom(1., 0.);
  complex<double> poly;
  vector<complex<double>> z(n);
  bool flag = false;
  while (flag == false) {
    flag = true;
    for (int k = 0; k < n; k++) {
      for (int j = 0; j < n; j++) {
        if (k != j) {
          denom *= (roots[k] - roots[j]);
        }
      }
      poly = Poly(n, sigma, roots[k]);
      z[k] = poly / denom;
      roots[k] = roots[k] - z[k];
      denom = 1.;
    }
    for (int m = 0; m < n; m++) {
      if (abs(z[m]) > accuracy) {
        flag = false;
      }
    }
  }
  return roots;
}
