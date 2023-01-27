#pragma once

#include <complex>
#include <fstream>
#include <iostream>

#include "accum_type.hh"
#include "dagger.hh"
#include "traceless_antiherm.hh"

using Complex = std::complex<double>;

class _u1 {
public:
  const size_t N_c = 1;
  explicit _u1() : a(0) {}
  explicit _u1(double _a) : a(_a) {}
  _u1(const _u1 &U) : a(U.a) {}
  _u1(Complex _a) : a(std::arg(_a)) {}

  friend Complex operator+(const _u1 &U1, const _u1 &U2);
  friend Complex operator-(const _u1 &U1, const _u1 &U2);
  friend _u1 operator*(const _u1 &U1, const _u1 &U2);

  // implicit conversion operator to complex
  operator Complex() const { return (std::exp(a * Complex(0., 1.))); }
  _u1 &operator*=(const _u1 &U1) {
    this->a += U1.a;
    return *this;
  }
  _u1 round(size_t n) const {
    double dn = n;
    return _u1(double(std::round(a) * dn) / dn);
  }

  double geta() const { return (a); }
  void operator=(const _u1 &U) { a = U.a; }
  void set(const double _a) { a = _a; }
  _u1 dagger() const { return (_u1(-a)); }
  double retrace() const { return (std::cos(a)); }
  Complex det() const { return (std::exp(a * Complex(0., 1.))); }
  void restoreSU() {}

private:
  double a;
};

inline double retrace(_u1 const &U) {
  return (std::cos(U.geta()));
}

inline double retrace(const Complex c) {
  return (std::real(c));
}

inline Complex trace(const Complex c) {
  return (c); // for U(1) the trace operator acts trivially
}

template <> struct accum_type<_u1> {
  typedef Complex type;
};

template <> inline Complex dagger(const Complex &x) {
  return std::conj(x);
}

template <> inline Complex traceless_antiherm(const Complex &x) {
  return (Complex(0., std::imag(x)));
}

//_u1 operator*(const _u1 &U1, const _u1 &U2);
// Complex operator*(const _u1 &U1, const Complex &U2);
// Complex operator*(const Complex &U1, const _u1 &U2);
// Complex operator+(const _u1 &U1, const _u1 &U2);
// Complex operator-(const _u1 &U1, const _u1 &U2);
// void operator+=(Complex & U1, const _u1 & U2);
// void operator*=(Complex & U1, const _u1 & U2);

inline _u1 operator*(const _u1 &U1, const _u1 &U2) {
  //  _u1 res;
  //  res.a = U1.a + U2.a;
  return _u1(U1.a + U2.a);
}

inline Complex operator*(const _u1 &U1, const Complex &U2) {
  return (std::exp(U1.geta() * Complex(0., 1.)) * U2);
}
inline Complex operator*(const Complex &U1, const _u1 &U2) {
  return (U1 * std::exp(U2.geta() * Complex(0., 1.)));
}

inline Complex operator+(const _u1 &U1, const _u1 &U2) {
  return (std::exp(U1.a * Complex(0., 1.)) + std::exp(U2.a * Complex(0., 1.)));
}

inline Complex operator-(const _u1 &U1, const _u1 &U2) {
  return (std::exp(U1.a * Complex(0., 1.)) - std::exp(U2.a * Complex(0., 1.)));
}

inline void operator+=(Complex &U1, const _u1 &U2) {
  U1 += std::exp(U2.geta() * Complex(0., 1.));
}

inline void operator*=(Complex &U1, const _u1 &U2) {
  U1 *= std::exp(U2.geta() * Complex(0., 1.));
}
