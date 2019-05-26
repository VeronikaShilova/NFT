#pragma once

#include "ProgramParams.h"
#include "GammaFunction.h"

// nonlinear spectrum of potential or signal
class Spectrum {
private:
  vector<complex<double>> disc_; // discrete component
  vector<complex<double>> cont_; // continuous component
  const double Xi1_ = -107.; // sufficienlty large interval
  const double Xi2_ = 107.; // boundaries for nonlinear signal frequency
public:
  double Xi1() const { return Xi1_; }
  double Xi2() const { return Xi2_; }
  void set_disc(vector<complex<double>> && v) {
    disc_ = forward<vector<complex<double>>>(v);
  }
  void set_cont(vector<complex<double>> && v) {
    cont_ = forward<vector<complex<double>>>(v);
  }
  const vector<complex<double>>& disc() const { return disc_; }
  const vector<complex<double>>& cont() const { return cont_; }
  double h() const {
    assert(cont_.size() != 0, "Continuous spectrum vector is empty, can't count h!");
    return (Xi2_ - Xi1_) / cont_.size();
  }
};

// potential or signal
class Potential { 
protected:
  const double A_; // amplitude of signal
  const double T1_; // signal is propagated over time span from T1 to T2
  const double T2_;
  Spectrum s_;
  Potential(const double A, const double T1, const double T2) : A_(A), T1_(T1), T2_(T2) {}
public:
  void set_disc_sp(vector<complex<double>> && v) {
    s_.set_disc(forward<vector<complex<double>>>(v));
  }
  void set_cont_sp(vector<complex<double>> && v) {
    s_.set_cont(forward<vector<complex<double>>>(v));
  }
  double A() const { return A_; }
  double T1() const { return T1_; }
  double T2() const { return T2_; }
  const vector<complex<double>>& disc_sp() const { return s_.disc(); }
  const vector<complex<double>>& cont_sp() const { return s_.cont(); }
  const double Xi1_sp() const { return s_.Xi1(); }
  const double Xi2_sp() const { return s_.Xi2(); }
  virtual complex<double> q(double t) const = 0;
  // a(xi), b(xi) are continuous spectral functions for a given signal
  virtual complex<double> a(const complex<double>& xi) const = 0;
  virtual complex<double> b(const complex<double>& xi) const = 0;
  virtual complex<double> r(const complex<double>& xi) const = 0;
  virtual complex<double> l(const complex<double>& xi) const = 0;
  // Fourrier transform of l(xi)=conj(b(xi))/a(xi)
  virtual complex<double> L(const double x) const = 0;
  virtual vector<complex<double>> get_analytical(const unsigned int N) const = 0;
};

// Satsuma-Yajima potential
class Satsuma : public Potential {
public:
  Satsuma(const double A, const double T1, const double T2) : Potential(A, T1, T2) {}
  complex<double> q(double t) const override {
    return complex<double>(A_ / cosh(t));
  }
  complex<double> a(const complex<double>& xi) const override {
    return complex<double>(pow(gamma(0.5 - i*xi, 0), 2) / (gamma(-i*xi + 0.5 + A_, 0)*gamma(0.5 - A_ - i*xi, 0)));
  }
  complex<double> b(const complex<double>& xi) const override {
    return complex<double>(-sin(pi*A_) / cosh(pi*xi));
  }
  complex<double> r(const complex<double>& xi) const override {
    return complex<double>(-sin(pi*A_)*gamma(-i*xi + 0.5 + A_, 0)*gamma(0.5 - A_ - i*xi, 0) / (pow(gamma(0.5 - i*xi, 0), 2)*cosh(pi*xi)));
  }
  complex<double> l(const complex<double>& xi) const {
    return complex<double>(conj(b(xi)) / a(xi));
  }
  complex<double> L(const double x) const override { // Fourrier transform of l(xi)
    complex<double> result = 0.;
    const double h = s_.h();
    vector<complex<double>> cont = s_.cont();
    for (int m = 0; m < cont.size(); m++) {
      result += cont[m] * exp(-i*(s_.Xi1() + (m + 0.5)*h)*x)*h;
    }
    result /= (2 * pi);
    vector<complex<double>> disc = s_.disc();
    for (int m = 0; m <= (disc.size()*0.5 - 1); m++) {
      result += -i*disc[disc.size()*0.5 + m] * exp(-i*disc[m] * x);
    }
    return result;
  }
  vector<complex<double>> get_analytical(const unsigned int N) const override {
    vector<complex<double>> result(N + 1);
    double h = (T2_ - T1_) / N;
    for (int j = 0; j <= N; j++) {
      result[j] = q(T1_ + j*h);
    }
    return result;
  }
};

// rectangular potential
class Rectangular : public Potential {
public:
  Rectangular(const double A, const double T1, const double T2) : Potential(A, T1, T2) {}
  complex<double> q(double t) const override {
    return complex<double>(A_);
  }
  complex<double> a(const complex<double>& xi) const override {
    complex<double> root = sqrt(pow(xi, 2) + pow(abs(A_), 2) + i);
    return complex<double>(exp(xi*i*(T2_ - T1_)))*(cos(root*(T2_ - T1_)) - (i*xi*(sin(root*(T2_ - T1_)))) / root);
  }
  complex<double> b(const complex<double>& xi) const override {
    complex<double> root = sqrt(pow(xi, 2) + pow(abs(A_), 2) + i);
    return complex<double>(-(A_ *sin(root*(T2_ - T1_))) / root);
  }
  complex<double> r(const complex<double>& xi) const override {
    return complex<double>(b(xi)/a(xi));
  }
  complex<double> l(const complex<double>& xi) const {
    return complex<double>(conj(b(xi)) / a(xi));
  }
  complex<double> L(const double x) const override {
    complex<double> result = 0.;
    const double h = s_.h();
    vector<complex<double>> cont = s_.cont();
    for (int m = 0; m < cont.size(); m++) {
      result += cont[m] * exp(-i * (s_.Xi1() + (m + 0.5) * h) * x) * h;
    }
    result /= (2 * pi);
    vector<complex<double>> disc = s_.disc();
    for (int m = 0; m <= (disc.size() * 0.5 - 1); m++) {
      result += -i * disc[disc.size() * 0.5 + m] * exp(-i * disc[m] * x);
    }
    return result;
  }
  vector<complex<double>> get_analytical(const unsigned int N) const override {
    vector<complex<double>> result(N + 1);
    double h = (T2_ - T1_) / N;
    for (int j = 0; j <= N; j++) {
      result[j] = q(T1_ + j*h);
    }
    return result;
  }
};

// multisoliton potential
class nSoliton : public Potential {
private:
  const int numberOfSolitions_;
  const int solitonLength_;
public:
  nSoliton(const double numberOfSolitions, const double solitonLength, const double solitonAmplitude) :
    Potential(solitonAmplitude, 0., solitonLength*numberOfSolitions), 
    numberOfSolitions_(numberOfSolitions),
    solitonLength_(solitonLength) {}
  complex<double> q(double t) const override {
    for (int k = 1; k <= numberOfSolitions_; k++) {
      if (t < k * solitonLength_) {
        return complex<double>(exp(-i * (t - k * solitonLength_ * 0.5)) / cosh(t - k * solitonLength_ * 0.5));
      }
    }
  }
  complex<double> a(const complex<double>& xi) const override {
    return complex<double>(0);
  }
  complex<double> b(const complex<double>& xi) const override {
    return complex<double>(0);
  }
  complex<double> r(const complex<double>& xi) const override {
    return complex<double>(0);
  }
  complex<double> l(const complex<double>& xi) const {
    return complex<double>(0);
  }
  complex<double> L(const double x) const override {
    complex<double> result = 0.;
    const double h = s_.h();
    vector<complex<double>> cont = s_.cont();
    for (int m = 0; m < cont.size(); m++) {
      result += cont[m] * exp(-i * (s_.Xi1() + (m + 0.5) * h) * x) * h;
    }
    result /= (2 * pi);
    vector<complex<double>> disc = s_.disc();
    for (int m = 0; m <= (disc.size() * 0.5 - 1); m++) {
      result += -i * disc[disc.size() * 0.5 + m] * exp(-i * disc[m] * x);
    }
    return result;
  }
  vector<complex<double>> get_analytical(const unsigned int N) const override {
    vector<complex<double>> result(N + 1);
    double h = (T2_ - T1_) / N;
    for (int j = 0; j <= N; j++) {
      result[j] = q(T1_ + j*h);
    }
    return result;
  }
};

// dark soliton signal
class darkSoliton : public Potential {
private:
  const double alpha_, nu_;
public:
  darkSoliton(const double alpha, const double nu, const double T1, const double T2)
    : Potential(sqrt(1. + nu*nu), T1, T2)
    , alpha_(alpha)
    , nu_(nu) {}
  complex<double> q(double t) const override {
    double num = -4.*alpha_*nu_*A_*(A_ - 1.);
    double den = (A_ - 1.)*(A_ - 1.)*exp(-2.*A_*alpha_* t) +
                 nu_*nu_*exp(2.*A_*alpha_*t);
    return complex<double>(num/den);
  }
  complex<double> a(const complex<double>& xi) const override {
    return complex<double>(xi + i * alpha_);
  }
  complex<double> b(const complex<double>& xi) const override {
    return complex<double>(conj(i * nu_ * alpha_));
  }
  complex<double> r(const complex<double>& xi) const override {
    return complex<double>(b(xi)/a(xi));
  }
  complex<double> l(const complex<double>& xi) const {
    return complex<double>(i * nu_ * alpha_ / (xi + i * alpha_));
  }
  complex<double> L(const double x) const override {
    return complex<double>(nu_*alpha_*exp(-alpha_*x));
  }
  vector<complex<double>> get_analytical(const unsigned int N) const override {
    vector<complex<double>> result(N + 1);
    double h = (T2_ - T1_) / N;
    for (int j = 0; j <= N; j++) {
      result[j] = q(T1_ + j*h);
    }
    return result;
  }
};