#pragma once

#include "ProgramParams.h"

// intial values for iteration process are uniformly distributed on a circle with radius (rootMax-epsilon)
inline complex<double> root(int j, int n, double R) { return R * exp(2 * pi * j / n); }

class Polynom {
private:
  vector<complex<double>> const& coefs_;
  double maxRoot_;
  inline complex<double> power(complex<double> z, unsigned int n) const {
    assert(n != 0);
    complex<double> res = z;
    for (int i = 1; i < n; i++)
      res *= z;
    return res;
  }
public:
  Polynom(vector<complex<double>> const& coefs) : coefs_(coefs) {
    maxRoot_ = abs(coefs_[0]);
    for (int k = 1; k < coefs_.size(); k++) {
      double a = abs(coefs_[k]);
      if (a > maxRoot_) maxRoot_ = a;
    }
    maxRoot_ += 1.;
  }
  inline complex<double> operator()(complex<double> z) const {
    size_t n = coefs_.size();
    complex<double> res = coefs_[n - 1];
    if (z != complex<double>(0., 0.)) {
      for (int k = 1; k < n; k++)
        res += coefs_[k] * power(z, n - k);
      res += power(z, n);
    }
    return res;
  }
  inline int    deg()     const { return coefs_.size(); }
  inline double rootMax() const { return maxRoot_; }
};

vector<complex<double>> findRoots(Polynom const& p) {
  int n = p.deg();
  double max = p.rootMax();
  const double eps = 1e-1;
  const double acc = 1e-10;
  auto r = [=](int j) { return root(j, n, max - eps); };
  vector<complex<double>> res(n);
  for (int k = 0; k < n; k++)
    res[k] = r(k);
  bool flag = false;
  while (flag == false) {
    flag = true;
    for (int k = 0; k < n; k++) {
      complex<double> denom(1., 0.);
      for (int j = 0; j < k; j++)
        denom *= (res[k] - res[j]);
      for (int j = k + 1; j < n; j++)
        denom *= (res[k] - res[j]);
      complex<double> zk = p(res[k]) / denom;
      if (abs(zk) > acc) flag = false;
      res[k] -= zk;
    }
  }
  return res;
}


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
