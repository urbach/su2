/**
 * @file flat_spacetime_gauge_energy.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief gauge energy in flat spacetime (euclidean metric)
 * @version 0.1
 * @date 2022-05-20
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include "gaugeconfig.hh"
#ifdef _USE_OMP_
#include <omp.h>
#endif

namespace flat_spacetime {

  /**
   * @brief real part of the trace of the Wilson plaquette
   *
   * \sum_mu \sum_nu<mu Re(Tr(P_{mu nu}))
   *
   * @tparam T
   * @param U
   * @param spatial_only: true when only the plaquettes with mu, nu > 0 are calculated
   * @param xi anisotropy in the action: S_G \supset (\beta/xi)*P_{0i} + (xi*\beta)*P_{ij}
   * @return double
   */
  template <class T>
  double retr_sum_Wplaquettes(const gaugeconfig<T> &U,
                              const double &xi = 1.0,
                              const bool &anisotropic = false,
                              const bool &spatial_only = false) {
    double res = 0.;
    size_t startmu = spatial_only; // 0 if spatial_only==false, 1 if spatial_only==true

    if (!anisotropic) {
#pragma omp parallel for reduction(+ : res)
      for (size_t x0 = 0; x0 < U.getLt(); x0++) {
        for (size_t x1 = 0; x1 < U.getLx(); x1++) {
          for (size_t x2 = 0; x2 < U.getLy(); x2++) {
            for (size_t x3 = 0; x3 < U.getLz(); x3++) {
              std::vector<size_t> x = {x0, x1, x2, x3};
              std::vector<size_t> xplusmu = x;
              std::vector<size_t> xplusnu = x;
              for (size_t mu = startmu; mu < U.getndims() - 1; mu++) {
                for (size_t nu = mu + 1; nu < U.getndims(); nu++) {
                  xplusmu[mu] += 1;
                  xplusnu[nu] += 1;
                  res += retrace(U(x, mu) * U(xplusmu, nu) * U(xplusnu, mu).dagger() *
                                 U(x, nu).dagger());
                  xplusmu[mu] -= 1;
                  xplusnu[nu] -= 1;
                }
              }
            }
          }
        }
      }
    }
    // 2n option - anisotropic lattice present
    if (anisotropic) {
#pragma omp parallel for reduction(+ : res)
      for (size_t x0 = 0; x0 < U.getLt(); x0++) {
        for (size_t x1 = 0; x1 < U.getLx(); x1++) {
          for (size_t x2 = 0; x2 < U.getLy(); x2++) {
            for (size_t x3 = 0; x3 < U.getLz(); x3++) {
              std::vector<size_t> x = {x0, x1, x2, x3};
              std::vector<size_t> xplusmu = x;
              std::vector<size_t> xplusnu = x;
              for (size_t mu = startmu; mu < U.getndims() - 1; mu++) {
                for (size_t nu = mu + 1; nu < U.getndims(); nu++) {
                  xplusmu[mu] += 1;
                  xplusnu[nu] += 1;
                  double eta = xi;
                  if ((mu==0) ^ (nu==0)) {
                    // at least one direction is temporal (but not both)
                    eta = 1.0 / xi;
                  }
                  res += eta * retrace(U(x, mu) * U(xplusmu, nu) *
                                       U(xplusnu, mu).dagger() * U(x, nu).dagger());
                  xplusmu[mu] -= 1;
                  xplusnu[nu] -= 1;
                }
              }
            }
          }
        }
      }
    }
    return res;
  }

  /** [DEPRECATED - USE retr_sum_Wplaquettes()]
   * @brief Wilson plaquette gauge energy
   *
   * \sum_mu \sum_nu<mu tr(P_{mu nu})
   *
   * checked for gauge invariance
   * if spatial_only, only the plaquettes with mu, nu > 0 are calculated
   *
   * @tparam T
   * @param U
   * @param spatial_only
   * @return double
   */
  template <class T>
  double gauge_energy(const gaugeconfig<T> &U, bool spatial_only = false) {
    double res = 0.;
    size_t startmu = spatial_only; // 0 if spatial_only==false, 1 if spatial_only==true

#pragma omp parallel for reduction(+ : res)
    for (size_t x0 = 0; x0 < U.getLt(); x0++) {
      for (size_t x1 = 0; x1 < U.getLx(); x1++) {
        for (size_t x2 = 0; x2 < U.getLy(); x2++) {
          for (size_t x3 = 0; x3 < U.getLz(); x3++) {
            std::vector<size_t> x = {x0, x1, x2, x3};
            std::vector<size_t> xplusmu = x;
            std::vector<size_t> xplusnu = x;
            for (size_t mu = startmu; mu < U.getndims() - 1; mu++) {
              for (size_t nu = mu + 1; nu < U.getndims(); nu++) {
                xplusmu[mu] += 1;
                xplusnu[nu] += 1;
                res += retrace(U(x, mu) * U(xplusmu, nu) * U(xplusnu, mu).dagger() *
                               U(x, nu).dagger());
                xplusmu[mu] -= 1;
                xplusnu[nu] -= 1;
              }
            }
          }
        }
      }
    }
    return res;
  }

} // namespace flat_spacetime