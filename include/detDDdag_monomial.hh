// detDDdag_monomial
/*
  Simone Romiti - simone.romiti@uni-bonn.de

  Evaluation of det(D*D^{\dagger}) with the pseudo-fermion action
  as in eq. (10) of https://www.sciencedirect.com/science/article/pii/0550321389903246
  or s eq. (8.9) of https://link.springer.com/book/10.1007/978-3-642-01850-3
*/

#pragma once

#include "adjointfield.hh"
#include "gauge_energy.hpp"
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "su2.hh"
#include "u1.hh"
#include <array>
#include <complex>
#include <random>
#include <vector>

#include "staggered.hpp" // spinor object

namespace staggered {
  template <class iT> using nd_max_arr = spacetime_lattice::nd_max_arr<iT>;

  // detDDdag monomial : evaluation of det(D*D^{\dagger}) through pseudo-fermions
  template <typename Float, class Group>
  class detDDdag_monomial : public monomial<Float, Group> {
    typedef std::complex<Float> Complex;

    size_t n_traj = 0;

  public:
    Float m0; // bare mass (in lattice units)

    std::string SOLVER; // type of the SOLVER
    Float TOLERANCE; // tolerance of the CG solver
    size_t VERBOSITY; // verbosity of the CG solver
    size_t SEED; // seed of the random number generator

    // pseudo-fermion field phi, kept constant along the MD trajectory
    staggered::spinor_lat<Float, Complex> phi;

    detDDdag_monomial<Float, Group>(unsigned int _timescale,
                                    const Float &m0_val,
                                    const std::string &solver,
                                    const Float &tolerance,
                                    const size_t &seed,
                                    const size_t &verb)
      : monomial<Float, Group>::monomial(_timescale) {
      m0 = m0_val;

      SOLVER = solver;
      TOLERANCE = tolerance;
      SEED = seed;
      VERBOSITY = verb;
    }

    // during the heatbath we generate R and store phi
    void heatbath(const hamiltonian_field<Float, Group> &h) override {
      const int N = h.U->getVolume(); // total number of lattice points

      const size_t Lt = h.U->getLt(), Lx = h.U->getLx(), Ly = h.U->getLy(),
                   Lz = h.U->getLz();
      const nd_max_arr<size_t> dims = {Lt, Lx, Ly, Lz}; // vactor of spacetime dimensions

      // generating the gaussian R

      const staggered::spinor_lat<Float, Complex> R =
        staggered::gaussian_spinor<Float, Complex>(
          dims, 0.0, 1.0 / sqrt(2), SEED + n_traj); // e^{-x^2} has sigma=1/sqrt(2)

      n_traj += 1;

      (*this).phi = staggered::apply_D<Float, Complex>(*h.U, (*this).m0, R);

      monomial<Float, Group>::Hold = R.norm_squared(); // R^{\dagger}*R is real
      return;
    }

    /**
     * @brief accept the new configuration and update the hamiltonian
     * After the MD trajectory we compute the new
     * R^{\dagger}*R = \phi[U_0]^{\dagger} * (D*D^{\dagger})[U_1]^{-1} * \phi[U_0]
     * where phi is the one computed at the beginning of the trajectory (configuration
     * U_0) and D, D^{\dagger} are evaluated from the new gauge config U_1
     * @param h hamiltonian field
     */
    void accept(const hamiltonian_field<Float, Group> &h) override {
      // Operator D*D^{\dagger} . Hermitian and invertible -> can apply the CG inversion
      const staggered::DDdag_matrix_lat<Float, Complex, _u1> DDdag(*h.U, (*this).m0);

      // applying (D*Ddag)^{-1} to \phi
      const staggered::spinor_lat<Float, Complex> chi =
        DDdag.inv((*this).phi, SOLVER, TOLERANCE, VERBOSITY, SEED);

      monomial<Float, Group>::Hnew = complex_dot_product(phi, chi).real();
      return;
    }

    /**
     * S_F  = \phi^{\dagger} * M^{-1} * \phi
     * dS_F = \phi^{\dagger} * M^{-1} * dM * M^{-1} * \phi = \chi^{\dagger} * dM * \chi
     * where \chi = M^{-1} * \phi
     * and M = D*D^{\dagger}
     * see eq. (8.44) of Gattringer&Lang
     */
    void derivative(adjointfield<Float, Group> &deriv,
                    const hamiltonian_field<Float, Group> &h,
                    const Float fac = 1.) const override {
      typedef typename accum_type<Group>::type accum;

      const staggered::DDdag_matrix_lat<Float, Complex, Group> DDdag(*h.U, (*this).m0);

      const staggered::spinor_lat<Float, Complex> chi =
        DDdag.inv((*this).phi, SOLVER, TOLERANCE, VERBOSITY, SEED);

      const size_t Lt = h.U->getLt(), Lx = h.U->getLx(), Ly = h.U->getLy(),
                   Lz = h.U->getLz();
      const size_t nd = h.U->getndims();

      const staggered::spinor_lat<Float, Complex> chi1 =
        staggered::apply_Ddag(*h.U, (*this).m0, chi);

      const Complex i(0.0, 1.0);

//#pragma omp target teams distribute parallel for //collapse(4)
// #pragma omp target teams distribute parallel for collapse(4)
#pragma omp parallel for // collapse(4)
      for (int x0 = 0; x0 < Lt; x0++) {
        for (int x1 = 0; x1 < Lx; x1++) {
          for (int x2 = 0; x2 < Ly; x2++) {
            for (int x3 = 0; x3 < Lz; x3++) {
              const nd_max_arr<int> x = {x0, x1, x2, x3};
              nd_max_arr<int> xm = x, xp = x;
              for (size_t mu = 0; mu < nd; mu++) {
                accum derSF;
                xm[mu]--; // x - mu
                xp[mu]++; // x + mu

                const Float eta_x_mu = staggered::eta(x, mu);
                const Float eta_xp_mu = staggered::eta(xp, mu);

                auto v_xp = (1.0 / 2.0) * eta_x_mu * (+i) * (*h.U)(x, mu) * chi1(xp);
                auto v_x =
                  -(1.0 / 2.0) * eta_xp_mu * (-i) * (*h.U)(x, mu).dagger() * chi1(x);

                // derivative of S_F with respect to U_{\mu}(x)
                derSF = conj(chi(x)) * v_xp + conj(chi(xp)) * v_x;
                derSF =
                  -2.0 * i * derSF; // get_deriv gives the imaginary part: Re(z) = Im(i*z)
                // std::cout << "derSF : " << derSF << "\n";
                deriv(x, mu) += get_deriv<Float>(derSF);

                xm[mu]++; // =x again
                xp[mu]--; // =x again
              }
            }
          }
        }
      }
      return;
    }
  };

  template <typename Float>
  class detDDdag_monomial<Float, _su2> : public monomial<Float, _su2> {
  public:
    detDDdag_monomial(unsigned int _timescale,
                      const Float &m0_val,
                      const std::string &solver,
                      const Float &tolerance,
                      const size_t &seed,
                      const size_t &verb)
      : monomial<Float, _su2>::monomial(_timescale) {
      spacetime_lattice::fatal_error("SU(2) not supported yet.", __func__);
    }

    void heatbath(const hamiltonian_field<Float, _su2> &h) override { return; }

    void accept(const hamiltonian_field<Float, _su2> &h) override { return; }

    void derivative(adjointfield<Float, _su2> &deriv,
                    const hamiltonian_field<Float, _su2> &h,
                    const Float fac = 1.) const override {
      return;
    }
  };

} // namespace staggered
