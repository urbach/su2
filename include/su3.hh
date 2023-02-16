/**
 * @file su3.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class and routines for an SU(3) matrix in the fundamental representation
 * @version 0.1
 * @date 2023-02-15
 *
 * @copyright Copyright (c) 2023
 *
 */

#pragma once

#include <complex>
#include <iostream>

#include "accum_type.hh"
#include "dagger.hh"
#include "traceless_antiherm.hh"

using Complex = std::complex<double>;

/**
 * @brief representation of an SU(2) matrix in the fundamental representation with the
 * convention of eq. (4.26) of Gattringer&Lang
 * https://link.springer.com/book/10.1007/978-3-642-01850-3
 *
 */
class _su3 {
public:
  const size_t N_c = 3;
  explicit _su3() {
    // default: identity matrix
    u = std::array<Complex, 3>({1.0, 0.0, 0.0});
    v = std::array<Complex, 3>({0.0, 1.0, 0.0});
  }
  explicit _su3(const std::array<Complex, 3> &_u, const std::array<Complex, 3> &_v)
    : u(_u), v(_v) {}
  _su3(const _su3 &U) : u(U.get_u()), v(U.get_v()) {}

  friend inline _su3 operator+(const _su3 &U1, const _su3 &U2);
  friend inline _su3 operator-(const _su3 &U1, const _su3 &U2);
  friend inline _su3 operator*(const _su3 &U1, const _su3 &U2);
  friend inline _su3 operator*(const Complex &U1, const _su3 &U2);
  friend inline _su3 operator*(const _su3 &U1, const Complex &U2);
  _su3 &operator*=(const _su3 &U2) {
    const std::array<Complex, 3> u1 = (*this).u;
    const std::array<Complex, 3> v1 = (*this).v;
    const std::array<Complex, 3> u2 = U2.u;
    const std::array<Complex, 3> v2 = U2.v;
    // vector components computed using sympy
    (*this).u = {u1[0] * u2[0] + u1[1] * v2[0] +
                   u1[2] * (std::conj(u2[1]) * std::conj(v2[2]) -
                            std::conj(u2[2]) * std::conj(v2[1])),
                 u1[0] * u2[1] + u1[1] * v2[1] +
                   u1[2] * (-std::conj(u2[0]) * std::conj(v2[2]) +
                            std::conj(u2[2]) * std::conj(v2[0])),
                 u1[0] * u2[2] + u1[1] * v2[2] +
                   u1[2] * (std::conj(u2[0]) * std::conj(v2[1]) -
                            std::conj(u2[1]) * std::conj(v2[0]))};
    (*this).v = {u2[0] * v1[0] + v1[1] * v2[0] +
                   v1[2] * (std::conj(u2[1]) * std::conj(v2[2]) -
                            std::conj(u2[2]) * std::conj(v2[1])),
                 u2[1] * v1[0] + v1[1] * v2[1] +
                   v1[2] * (-std::conj(u2[0]) * std::conj(v2[2]) +
                            std::conj(u2[2]) * std::conj(v2[0])),
                 u2[2] * v1[0] + v1[1] * v2[2] +
                   v1[2] * (std::conj(u2[0]) * std::conj(v2[1]) -
                            std::conj(u2[1]) * std::conj(v2[0]))};
    return *this;
  }
  _su3 round(size_t n) const {
    double dn = n;
    std::array<Complex, 3> u2, v2;
    for (size_t i = 0; i < N_c; i++) {
      u2[i] = Complex(std::round(std::real(u[i]) * dn) / dn,
                      std::round(std::imag(u[i]) * dn) / dn);
      v2[i] = Complex(std::round(std::real(v[i]) * dn) / dn,
                      std::round(std::imag(v[i]) * dn) / dn);
    }
    return _su3(u2, v2);
  }

  inline std::array<Complex, 3> get_u() const { return u; }
  inline std::array<Complex, 3> get_v() const { return v; }
  inline void operator=(const _su3 &U) {
    u = U.get_u();
    v = U.get_v();
  }
  inline void operator+=(const _su3 &U) {
    for (size_t i = 0; i < N_c; i++) {
      u[i] += U.u[i];
      v[i] += U.v[i];
    }
  }
  inline void operator-=(const _su3 &U) {
    for (size_t i = 0; i < N_c; i++) {
      u[i] -= U.u[i];
      v[i] -= U.v[i];
    }
  }
  void set(const std::array<Complex, 3> &_u, const std::array<Complex, 3> &_v) {
    u = _u;
    v = _v;
  }
  inline _su3 dagger() const {
    // vectors computed using sympy
    const std::array<Complex, 3> u2 = {std::conj(u[0]), std::conj(v[0]),
                                       u[1] * v[2] - u[2] * v[1]};
    const std::array<Complex, 3> v2 = {std::conj(u[1]), std::conj(v[1]),
                                       -u[0] * v[2] + u[2] * v[0]};

    return _su3(u2, v2);
  }
  inline Complex trace() const {
    // expression computed using sympy
    const Complex res =
      u[0] + v[1] + std::conj(u[0]) * std::conj(v[1]) - std::conj(u[1]) * std::conj(v[0]);
    return res;
  }
  inline double retrace() const { return std::real(this->trace()); }
  Complex det() const {
    // expression computed using sympy
    Complex res = u[0] * v[1] * std::conj(u[0]) * std::conj(v[1]) -
                  u[0] * v[1] * std::conj(u[1]) * std::conj(v[0]) +
                  u[0] * v[2] * std::conj(u[0]) * std::conj(v[2]) -
                  u[0] * v[2] * std::conj(u[2]) * std::conj(v[0]) -
                  u[1] * v[0] * std::conj(u[0]) * std::conj(v[1]) +
                  u[1] * v[0] * std::conj(u[1]) * std::conj(v[0]) +
                  u[1] * v[2] * std::conj(u[1]) * std::conj(v[2]) -
                  u[1] * v[2] * std::conj(u[2]) * std::conj(v[1]) -
                  u[2] * v[0] * std::conj(u[0]) * std::conj(v[2]) +
                  u[2] * v[0] * std::conj(u[2]) * std::conj(v[0]) -
                  u[2] * v[1] * std::conj(u[1]) * std::conj(v[2]) +
                  u[2] * v[1] * std::conj(u[2]) * std::conj(v[1]);

    return res;
  }
  /**
   * @brief formula (4.27) of Gattringer&Lang
   * https://link.springer.com/book/10.1007/978-3-642-01850-3
   *
   */
  void restoreSU() {
    double n_u = 0.0; // norms of 'u' and 'v2'
    for (size_t i = 0; i < N_c; i++) {
      n_u += std::pow(std::abs((*this).u[i]), 2.0);
    }
    n_u = std::sqrt(n_u);
    for (size_t i = 0; i < N_c; i++) {
      (*this).u[i] /= n_u;
    }
    const Complex v_ustar =
      v[0] * std::conj(u[0]) + v[1] * std::conj(u[1]) + v[2] * std::conj(u[2]);
    std::array<Complex, 3> v2 = (*this).v;
    double n_v2 = 0.0;
    for (size_t i = 0; i < N_c; i++) {
      v2[i] -= v_ustar * u[i];
    }
    for (size_t i = 0; i < N_c; i++) {
      n_v2 += std::pow(std::abs(v2[i]), 2.0);
    }
    n_v2 = std::sqrt(n_v2);
    for (size_t i = 0; i < N_c; i++) {
      v2[i] /= n_v2;
    }
    (*this).v = v2;
    return;
  }

private:
  std::array<Complex, 3> u, v;
};

inline Complex trace(_su3 const &U) {
  return U.trace();
}

inline double retrace(_su3 const &U) {
  return U.retrace();
}

template <> inline _su3 dagger(const _su3 &u) {
  return u.dagger();
}

inline _su3 operator*(const _su3 &U1, const _su3 &U2) {
  _su3 U = U1;
  U *= U2;
  return U;
}

inline _su3 operator+(const _su3 &U1, const _su3 &U2) {
  _su3 U = U1;
  U += U2;
  return U;
}

inline _su3 operator-(const _su3 &U1, const _su3 &U2) {
  _su3 U = U1;
  U -= U2;
  return U;
}

inline _su3 operator*(const Complex &z, const _su3 &U2) {
  const size_t N_c = U2.N_c;
  std::array<Complex, 3> u, v;
  for (size_t i = 0; i < N_c; i++) {
    u[i] *= z;
    v[i] *= z;
  }
  return _su3(u, v);
}

inline _su3 operator*(const _su3 &U1, const Complex &z) {
  return z * U1;
}

template <> inline _su3 traceless_antiherm(const _su3 &x0) {
  const _su3 Id({1.0, 0.0, 0.0},
                {0.0, 1.0, 0.0}); // default constructor: identity operator
  _su3 x = x0;
  x = x - (trace(x) / double(x.N_c)) * Id;
  x = 0.5 * (x - x.dagger());
  return x;
}

using su3 = _su3;