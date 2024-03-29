#pragma once

#include "gaugeconfig.hh"

/**
 * @brief Update the value of the staple attached to U_{\mu}(x)
 * Note: use this function only once after the staple object has been created.
 * This function updates the staples surrounding the link at position x and direction mu:
 * U^staples = \sum_{nu != mu} U_nu(x+mu) U_mu(x+nu)^dagger + U_nu(x)^dagger
 *             + sum__{nu != mu} U_nu(x+mu-nu)^dagger U_mu(x-nu)^dagger U_nu(x-nu)
 * For anisotropic lattices, the action requires weighting the temporal and spacial
 * plaquettes differently, see e.g. https://arxiv.org/abs/hep-lat/0209159. This is
 * achieved by including this factor in the staples For other cases, only the staples
 * involving only spatial links are needed, then the summation starts at nu=1 to exclude
 * temporal links. This is done by setting spatial_only to true
 * @tparam T
 * @tparam S symmetry group
 * @param K staple
 * @param U gauge config
 * @param x spacetime point of
 * @param mu spacetime direction
 * @param xi anisotropy parameter
 * @param anisotropic boolean flag: true when considering anisotropy
 * @param spatial_only boolean flag: true when summing over spatial staples only
 */
template <class T, class S, class Arr>
void get_staples(T &K,
                 gaugeconfig<S> &U,
                 Arr const x,
                 const size_t mu,
                 const double xi = 1.0,
                 bool anisotropic = false,
                 bool spatial_only = false) {
  size_t startnu = 0;
  if (spatial_only) {
    startnu = 1;
  }
  Arr x1 = x, x2 = x;
  x1[mu] += 1;
  if (!anisotropic) {
    for (size_t nu = startnu; nu < U.getndims(); nu++) {
      if (nu != mu) {
        x2[nu]++;
        K += U(x1, nu) * U(x2, mu).dagger() * U(x, nu).dagger();
        x2[nu]--;
      }
    }
    for (size_t nu = startnu; nu < U.getndims(); nu++) {
      if (nu != mu) {
        x1[nu]--;
        x2[nu]--;
        K += U(x1, nu).dagger() * U(x2, mu).dagger() * U(x2, nu);
        x2[nu]++;
        x1[nu]++;
      }
    }
  }
  if (anisotropic) {
    double factor;
    for (size_t nu = startnu; nu < U.getndims(); nu++) {
      if (nu != mu) {
        factor = (((nu == 0) || (mu == 0)) ? 1.0 / xi : xi);
        x2[nu]++;
        K += factor * U(x1, nu) * U(x2, mu).dagger() * U(x, nu).dagger();
        x2[nu]--;
      }
    }
    // Maybe put both loops together so only one ?: operator is needed?
    for (size_t nu = startnu; nu < U.getndims(); nu++) {
      if (nu != mu) {
        factor = (((nu == 0) || (mu == 0)) ? 1.0 / xi : xi);
        x1[nu]--;
        x2[nu]--;
        K += factor * U(x1, nu).dagger() * U(x2, mu).dagger() * U(x2, nu);
        x2[nu]++;
        x1[nu]++;
      }
    }
  }
}

/**
 * @brief Get the staples for the APE smearing
 * see eq. (13) of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.70.014504
 * staples have to be 'daggered' / in the other direction compared to staples for action:
 * action forms plaquette, smearing needs staples parallel to link
 */
template <class T, class S, class Arr>
void get_staples_APE(
  T &K, const gaugeconfig<S> &U, Arr const x, const size_t mu, bool spatial_only = false) {
  size_t startnu = size_t(spatial_only); // 0 or 1, casted from bool
  Arr y = x, z = x;
  for (size_t nu = startnu; nu < U.getndims(); nu++) {
    if (nu != mu) {
      y[mu]++; // y = x+\mu
      z[nu]++; // z = x + \nu
      K += U(x, nu) * U(z, mu) * U(y, nu).dagger();

      z[nu] -= 2; // z = x - \nu
      y[mu]--; // y = x

      y[nu]--; // y = x - \nu
      z[mu]++; // z = x + \mu -\nu
      K += U(y, nu).dagger() * U(y, mu) * U(z, nu);
      y[nu]++; // y = x
      z[mu]--; // z = x - \nu
      z[nu]++; // z = x
    }
  }
}
