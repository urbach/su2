/**
 * @file gaugemonomial.hh
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief gauge monomial S_G - piece of the action containing gluons only
 * @version 0.1
 * @date 2022-05-26
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once
#include "adjointfield.hh"
#include "gauge_energy.hpp"
#include "gaugeconfig.hh"
#include "geometry.hh"
#include "get_staples.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "su2.hh"
#include "u1.hh"

#include <complex>
#include <vector>

namespace flat_spacetime {

  /**
   * @brief Get the S_G action for the HMC
   * Get the action subtracted from unnecessary constant factors (unphisical).
   *
   * Returns:
   * (\beta/Nc)*\sum_{x} \sum_{\mu<\nu} \eta_{\mu \nu} Re(Tr(P_{\mu \nu}(x)))
   * with \eta_{0i}=1/\xi and \eta_{ij}=\xi
   *
   * where beta = 2*N_c/g_0^2
   *
   * @tparam Float
   * @tparam Group
   * @param U gauge configuration
   * @param xi anisotropy (==1.0 when the lattice is isotropic)
   * @param anisotropic true when considering the anisotropy
   * @return Float
   */
  template <typename Float, class Group>
  Float get_S_G_hmc(const gaugeconfig<Group> &U,
                    const double &xi = 1.0,
                    const bool &anisotropic = false) {
    const Float res = -U.getBeta() *
                      flat_spacetime::retr_sum_Wplaquettes(U, xi, anisotropic) /
                      double(U.getNc());
    return res;
  }

  /**
   * @brief Full gauge action
   * Returns
   *
   * S_G = (\beta/Nc)*\sum_{x} \sum_{\mu<\nu}
   * \eta_{\mu \nu} Re(Tr(1 - P_{\mu\nu}(x)))
   *
   * where beta = 2*N_c/g_0^2
   *
   * @tparam Float
   * @tparam Group
   * @param U
   * @param xi
   * @param anisotropic
   * @return Float
   */
  template <typename Float, class Group>
  Float get_S_G(const gaugeconfig<Group> &U,
                const double &xi = 1.0,
                const bool &anisotropic = false) {
    const size_t Dims_fact = spacetime_lattice::num_pLloops_half(U.getndims());
    const Float res = U.getBeta() * (U.getVolume() * Dims_fact) +
                      flat_spacetime::get_S_G_hmc<Group>(U, xi, anisotropic);
    return res;
  }

  // gauge monomial
  template <typename Float, class Group>
  class gaugemonomial : public monomial<Float, Group> {
  public:
    gaugemonomial<Float, Group>(unsigned int _timescale, const double &_xi = 1.0)
      : monomial<Float, Group>::monomial(_timescale) {
      xi = _xi;
      if (_xi != 1.0) {
        anisotropic = true;
      }
    }
    // S_g = sum_x sum_{mu<nu} beta*(1- 1/Nc*Re[Tr[U_{mu nu}]])
    // beta = 2*N_c/g_0^2
    void heatbath(hamiltonian_field<Float, Group> const &h) override {
      monomial<Float, Group>::Hold =
        flat_spacetime::get_S_G_hmc<Float, Group>(*h.U, (*this).xi, (*this).anisotropic);
      return;
    }
    void accept(hamiltonian_field<Float, Group> const &h) override {
      monomial<Float, Group>::Hnew =
        flat_spacetime::get_S_G_hmc<Float, Group>(*h.U, (*this).xi, (*this).anisotropic);
      return;
    }

    /**
     * @brief updating the derivative of the action with respect to the gauge field
     * contribution in this monomial
     *
     * @param deriv derivative to be updated
     * @param h hamiltonian field (Hybrid Monte Carlo evolution)
     * @param fac multiplicative facator of the derivative (1.0 for the HMC but != 1 in
     * other cases e.g. the gradient flow evolution)
     */
    void derivative(adjointfield<Float, Group> &deriv,
                    hamiltonian_field<Float, Group> const &h,
                    const Float fac = 1.) const override {
      typedef typename accum_type<Group>::type accum;
#pragma omp parallel for
      for (size_t x0 = 0; x0 < h.U->getLt(); x0++) {
        for (size_t x1 = 0; x1 < h.U->getLx(); x1++) {
          for (size_t x2 = 0; x2 < h.U->getLy(); x2++) {
            for (size_t x3 = 0; x3 < h.U->getLz(); x3++) {
              std::vector<size_t> x = {x0, x1, x2, x3};
              for (size_t mu = 0; mu < h.U->getndims(); mu++) {
                accum S;
                get_staples(S, *h.U, x, mu, (*this).xi, (*this).anisotropic);
                S = (*h.U)(x, mu) * S; // U*A in eq. 8.40 in Gattringer&Lang

                const double num_fact_i = fac * h.U->getBeta() / double(h.U->getNc());
                deriv(x, mu) += num_fact_i * get_deriv<double>(S);
              }
            }
          }
        }
      }
      return;
    }

  private:
    // size_t Dims_fact; // d*(d-1)/2

    bool anisotropic = false;
    double xi; // bare anisotropy xi
  };

} // namespace flat_spacetime
