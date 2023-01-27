// staggered.hpp
/**
 Simone Romiti - simone.romiti@uni-bonn.de

 Staggered fermions routines

 Note: Instead of std::vector it is used nd_max_arr<int>, because we always deal
 with 4-dimensional containes for coordinates and dimensions
*/

#pragma once

#include <array>
#include <complex>
#include <random>
#include <vector>

#include "adjointfield.hh"
#include "gauge_energy.hpp"
#include "gaugeconfig.hh"
#include "geometry.hh"
#include "get_staples.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "su2.hh"
#include "u1.hh"

#include "solver_type.hh" // generic type of solver

// given A(matrix) and b(vector), find x = A^{-1}*b
#include "BiCGStab.hpp" // Bi Conjugate Gradient
#include "CG.hpp" // standard conjugate gradient

namespace staggered {

  const size_t nd_max = spacetime_lattice::nd_max;
  template <class iT> using nd_max_arr = spacetime_lattice::nd_max_arr<iT>;

  // eta_{\mu}(x) as in eq. (16) of https://arxiv.org/pdf/2112.14640.pdf
  // or eq. (8) of https://www.sciencedirect.com/science/article/pii/0550321389903246
  double eta(const nd_max_arr<int> &x, const size_t &mu) {
    double s = 0;
    // TO BE TESTED #pragma omp parallel for reduction(+: s)
    for (int nu = 0; nu < mu; nu++) {
      s += x[nu];
    }
    return std::pow(-1.0, s);
  }

  // index of the lattice point given the dimensions
  // #pragma omp declare target
  size_t txyz_to_index(const size_t &t,
                       const size_t &x,
                       const size_t &y,
                       const size_t &z,
                       const size_t &Lt,
                       const size_t &Lx,
                       const size_t &Ly,
                       const size_t &Lz) {
    const geometry g(Lx, Ly, Lz, Lt); // note the order
    return g.getIndex(t, x, y, z);
  }
  // #pragma omp end declare target

  // overload : x={x0,x1,x2,x3}, dims={Lt,Lx,Ly,Lz}
  // #pragma omp declare target
  size_t txyz_to_index(const nd_max_arr<int> &x, const nd_max_arr<size_t> &dims) {
    return txyz_to_index(x[0], x[1], x[2], x[3], dims[0], dims[1], dims[2], dims[3]);
  }
  // #pragma omp end declare target

  /**
   * @brief vector of staggered "spinors" (no Dirac structure) for all the points of the
   * lattice Example: spinor_lat<double, std::complex<double> psi;
   * @tparam Float
   * @tparam Type
   */
  template <class Float, class Type> class spinor_lat {
  private:
    std::vector<Type> Psi;
    nd_max_arr<size_t> dims; // spacetime dimensions : {Lt, Lx, Ly, Lz}

  public:
    spinor_lat() {}
    ~spinor_lat() {}

    Type &operator()(const size_t &i) { return Psi[i]; }
    Type operator()(const size_t &i) const { return Psi[i]; }

    Type &operator()(const nd_max_arr<int> &x) {
      const size_t i = txyz_to_index(x, dims);
      return Psi[i];
    }

    Type operator()(const nd_max_arr<int> &x) const {
      const size_t i = txyz_to_index(x, dims);
      return Psi[i];
    }

    Type &operator[](const size_t &i) { return Psi[i]; }
    Type operator[](const size_t &i) const { return Psi[i]; }

    spinor_lat(const nd_max_arr<size_t> &_dims) {
      dims = _dims;
      const size_t n = spacetime_lattice::Npts_from_dims(dims);
      Psi.resize(n);
    }

    spinor_lat(const nd_max_arr<size_t> &_dims, const Type &val) {
      dims = _dims;
      const size_t n = spacetime_lattice::Npts_from_dims(dims);
      std::vector<Type> v(n, val);
      (*this).Psi = v;
    }

    nd_max_arr<size_t> get_dims() const { return dims; }

    size_t size() const { return Psi.size(); }

    spinor_lat<Float, Type> operator/(const Type &lambda) {
      const int N = this->size();
      spinor_lat<Float, Type> phi(this->get_dims());
      for (size_t i = 0; i < N; i++) {
        phi[i] = Psi[i] / lambda;
      }
      return phi;
    }

    void operator+=(const spinor_lat<Float, Type> &psi) { (*this) = (*this) + psi; }

    Type dot(const spinor_lat &b) const {
      return complex_dot_product((*this), b); // (*this)^{\dagger} * b
    }

    Float norm_squared() const {
      return complex_dot_product((*this), (*this)).real(); // A^{\dagger} * A is real
    }

    Float norm() const { return sqrt(this->norm_squared()); }

  }; // class spinor_lat

  // a + lambda*b
  template <class Float, class Type, class Type_lambda>
  spinor_lat<Float, Type> a_plus_lambda_b(const spinor_lat<Float, Type> &a,
                                          const Type_lambda &lambda,
                                          const spinor_lat<Float, Type> &b) {
    const int N = a.size();
    spinor_lat<Float, Type> c(a.get_dims());
    for (size_t i = 0; i < N; i++) {
      c[i] = a[i] + lambda * b[i];
    }
    return c;
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> operator+(const spinor_lat<Float, Type> &a,
                                    const spinor_lat<Float, Type> &b) {
    return a_plus_lambda_b(a, 1.0, b);
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> operator-(const spinor_lat<Float, Type> &a,
                                    const spinor_lat<Float, Type> &b) {
    return a_plus_lambda_b(a, -1.0, b);
  }

  // change of sign : psi --> -psi
  template <class Float, class Type>
  spinor_lat<Float, Type> operator-(const spinor_lat<Float, Type> &psi) {
    const spinor_lat<Float, Type> v(psi.get_dims());
    return (v - psi);
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> gaussian_spinor(const nd_max_arr<size_t> &dims,
                                          const Float &avr,
                                          const Float &sigma,
                                          const size_t &seed) {
    std::normal_distribution<Float> dis{avr, sigma};
    std::mt19937 gen(seed); //, gen_im(seed + 1);

    // Note: 'dims' is not used directly here, but is part of the correct initialization
    // of the spinor
    spinor_lat<Float, Type> psi_gauss(dims);
    const size_t n = psi_gauss.size();

    for (size_t i = 0; i < n; i++) { // lattice points
      const Type x = dis(gen);
      psi_gauss[i] = x;
    }
    return psi_gauss;
  }

  // \sum_{i} A_i^{\dagger}*B_i
  template <class Float, class Type>
  Type complex_dot_product(const spinor_lat<Float, Type> &A,
                           const spinor_lat<Float, Type> &B) {
    const int N = A.size();
    Type sum = 0.0;

    // #pragma omp parallel for
    for (size_t i = 0; i < N; i++) {
      sum += conj(A[i]) * B[i];
    }

    return sum;
  }

  template <class Float, class Type>
  Type operator*(const spinor_lat<Float, Type> &A, const spinor_lat<Float, Type> &B) {
    return complex_dot_product(A, B);
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> operator*(const Type &lambda,
                                    const spinor_lat<Float, Type> &psi) {
    const spinor_lat<Float, Type> v(psi.get_dims());
    return a_plus_lambda_b(v, lambda, psi);
  }

  template <class Float, class Type> Float norm(const spinor_lat<Float, Type> &psi) {
    return psi.norm();
  }

  // class for D^{\dagger}*D
  template <class Float, class Type, class Group> class DDdag_matrix_lat {
  public:
    gaugeconfig<Group> U;
    Float m;

    DDdag_matrix_lat() {}
    ~DDdag_matrix_lat() {}

    DDdag_matrix_lat(const gaugeconfig<Group> _U, const Float &_m) {
      U = _U;
      m = _m;
    }

    void operator=(const DDdag_matrix_lat &dd) {
      U = dd.U;
      m = dd.m;
    }

    // rows = cols = N = number of lattice points
    size_t rows() const { return U.getVolume(); }
    size_t cols() const { return this->rows(); }

    spinor_lat<Float, Type> inv(const spinor_lat<Float, Type> &psi,
                                const std::string &solver,
                                const Float &tol,
                                const size_t &verb,
                                const size_t &seed) const {
      typedef spinor_lat<Float, Type> LAvector;

      typedef DDdag_matrix_lat<Float, Type, Group> LAmatrix;

      typedef CG::LinearCG<Float, Type, LAmatrix, LAvector> cg;
      typedef typename solver_type<cg>::type svr_type;
      if (solver == "BiCGStab") {
        typedef BiCGStab::LinearBiCGStab<Float, Type, LAmatrix, LAvector> bcgstab;
        typedef typename solver_type<bcgstab>::type svr_type;
      }

      svr_type SVR((*this), psi);

      const LAvector phi0 =
        staggered::gaussian_spinor<Float, Complex>(psi.get_dims(), 0.0, 10.0, seed);

      if (verb > 1) {
        std::cout << "Calling the " << solver << " solver.\n";
      }

      SVR.solve(phi0, tol, verb);
      return SVR.get_solution();
    }
  };

  /**
   * @brief (D*D^{\dagger})*psi
   *
   * @tparam Float : double
   * @tparam Type : std::complex<double>
   * @tparam Group : su1 or su2
   * @param M : wrapper for D*D^{\dagger} matrix.
   *  Contains the information only about the configuration U and the quark mass
   * @param psi : spinor to which we apply D*D^{\dagger}
   * @return spinor_lat<Float, Type>
   */
  template <class Float, class Type, class Group>
  spinor_lat<Float, Type> operator*(const DDdag_matrix_lat<Float, Type, Group> &M,
                                    const spinor_lat<Float, Type> &psi) {
    gaugeconfig<Group> U = M.U;
    const Float m = M.m;

    const size_t Lt = U.getLt(), Lx = U.getLx(), Ly = U.getLy(), Lz = U.getLz();
    const size_t nd = U.getndims();
    const nd_max_arr<size_t> dims = psi.get_dims(); // vector of dimensions

    const int N = psi.size();
    return apply_D(U, m, apply_Ddag(U, m, psi));
  }

  // returns the result of D*psi, where D is the Dirac operator
  template <class Float, class Type, class _u1>
  spinor_lat<Float, Type> apply_D(const gaugeconfig<_u1> &U,
                                  const Float &m,
                                  const spinor_lat<Float, Type> &psi) {
    const size_t Lt = U.getLt(), Lx = U.getLx(), Ly = U.getLy(), Lz = U.getLz();
    const size_t nd = U.getndims();
    const nd_max_arr<size_t> dims = psi.get_dims(); // vector of dimensions

    const int N = psi.size();
    spinor_lat<Float, Type> phi(dims);

//#pragma omp target teams distribute parallel for //collapse(4)
#pragma omp parallel for
    for (int x0 = 0; x0 < Lt; x0++) {
      for (int x1 = 0; x1 < Lx; x1++) {
        for (int x2 = 0; x2 < Ly; x2++) {
          for (int x3 = 0; x3 < Lz; x3++) {
            const nd_max_arr<int> x = {x0, x1, x2, x3};
            nd_max_arr<int> xm = x, xp = x;
            for (size_t mu = 0; mu < nd; mu++) {
              const Float eta_x_mu = eta(x, mu);
              xm[mu]--; // x - mu
              xp[mu]++; // x + mu

              phi(x) += +(1.0 / 2.0) * eta_x_mu * U(x, mu) * psi(xp);
              phi(x) += -(1.0 / 2.0) * eta_x_mu * U(xm, mu).dagger() * psi(xm);

              xm[mu]++; // =x again
              xp[mu]--; // =x again
            }
            phi(x) += m * psi(x);
          }
        }
      }
    }
    return phi;
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> apply_D(const gaugeconfig<_su2> &U,
                                  const Float &m,
                                  const spinor_lat<Float, Type> &psi) {
    spacetime_lattice::fatal_error("Staggered fermions not supported for SU(2)",
                                   __func__);
    return psi;
  }

  // returns the reslut of D^{\dagger}*psi, where D is the Dirac operator
  template <class Float, class Type, class _u1>
  spinor_lat<Float, Type> apply_Ddag(const gaugeconfig<_u1> &U,
                                     const Float &m,
                                     const spinor_lat<Float, Type> &psi) {
    const size_t Lt = U.getLt(), Lx = U.getLx(), Ly = U.getLy(), Lz = U.getLz();
    const size_t nd = U.getndims();
    const nd_max_arr<size_t> dims = psi.get_dims(); // vector of dimensions

    const int N = psi.size();
    spinor_lat<Float, Type> phi(dims);
#pragma omp parallel for
    for (int x0 = 0; x0 < Lt; x0++) {
      for (int x1 = 0; x1 < Lx; x1++) {
        for (int x2 = 0; x2 < Ly; x2++) {
          for (int x3 = 0; x3 < Lz; x3++) {
            const nd_max_arr<int> x = {x0, x1, x2, x3};
            nd_max_arr<int> xm = x, xp = x;
            for (size_t mu = 0; mu < nd; mu++) {
              xm[mu]--; // x - mu
              xp[mu]++; // x + mu

              const Float eta_xm_mu = eta(xm, mu);
              const Float eta_xp_mu = eta(xp, mu);

              phi(x) += +(1.0 / 2.0) * eta_xm_mu * U(xm, mu).dagger() * psi(xm);
              phi(x) += -(1.0 / 2.0) * eta_xp_mu * U(x, mu) * psi(xp);

              xm[mu]++; // =x again
              xp[mu]--; // =x again
            }
            phi(x) += m * psi(x);
          }
        }
      }
    }
    return phi;
  }

  template <class Float, class Type>
  spinor_lat<Float, Type> apply_Ddag(const gaugeconfig<_su2> &U,
                                     const Float &m,
                                     const spinor_lat<Float, Type> &psi) {
    spacetime_lattice::fatal_error("Staggered fermions not supported for SU(2)",
                                   __func__);
    return psi;
  }

} // namespace staggered
