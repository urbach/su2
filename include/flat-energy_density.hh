#pragma once
#include "accum_type.hh"
#include "gaugeconfig.hh"
#include "tensors.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace flat_spacetime {

  /**
   * @brief symmetric definition of the energy density using the clover leaf
   *
   *
   *  <--   <--
   * |  ^  |  ^
   * |  |  |  | mu
   * -->   -->
   *
   *  <--   <--
   * |  ^  |  ^
   * |  |  |  |
   * -->   -->
   *        nu
   *
   * if !cloverdef only a single plaquette ist used instead of the
   * sum over the clover leafs
   *
   * checked for gauge invariance!
   *
   * energy density
   * E = 1/4 G_{mu nu}^a G_{mu nu}^a = 1/2 tr(G_{mu nu} G_{mu nu})
   *
   * topological charge
   * Q = 1./(32 pi^2) eps_\mu\nu\rho\sigma Trace[ G_\mu\nu G_\rho\sigma]
   *
   * from hep-lat/9603008 we take equation (6)
   * G_\mu\nu = 1/4 sum_clover 1/2 (clover - h.c.)_\mu\nu
   *
   * that means
   * Q = 1./(32 pi^2)/16 eps_\mu\nu\rho\sigma Trace[
   *  sum_clover 1/2 (clover - h.c.)_{\mu\nu} * sum_clover 1/2 (clover -
   * h.c.)_{\rho\sigma}]
   *
   * if we take only the terms with \mu < \nu and \rho < \sigma, we need
   * to multiply by a factor of 4. All 4 terms come with the same sign.
   *
   * @tparam T
   * @param U
   * @param res
   * @param Q
   * @param cloverdef
   * @param ss compute only the spatial-spatial contribution (continuum limit given by
   * eq. 3.1 of https://arxiv.org/pdf/1205.0781.pdf)
   */
  template <class T>
  void energy_density(const gaugeconfig<T> &U,
                      double &res,
                      double &Q,
                      bool cloverdef = true,
                      const bool &ss = false) {
    const size_t mu_start = bool(ss);

    res = 0.;
    Q = 0.;

    typedef typename accum_type<T>::type accum;
    // Euclidean 4D totally anti-symemtric tensor
    static epsilon4_t eps4 = new_epsilon4();

    std::vector<size_t> x = {0, 0, 0, 0};
    for (x[0] = 0; x[0] < U.getLt(); x[0]++) {
      for (x[1] = 0; x[1] < U.getLx(); x[1]++) {
        for (x[2] = 0; x[2] < U.getLy(); x[2]++) {
          for (x[3] = 0; x[3] < U.getLz(); x[3]++) {
            std::vector<size_t> x1 = x;
            std::vector<size_t> x2 = x;
            std::vector<size_t> x3 = x;
            accum G[4][4];
            for (size_t mu = mu_start; mu < U.getndims() - 1; mu++) {
              for (size_t nu = mu + 1; nu < U.getndims(); nu++) {
                x1[mu] += 1; // x + mu
                x2[nu] += 1; // x + nu
                accum leaf =
                  U(x, mu) * U(x1, nu) * U(x2, mu).dagger() * U(x, nu).dagger();
                x1[mu] -= 1; // x
                x2[nu] -= 1; // x

                if (cloverdef) {
                  x1[mu] -= 1; // x - mu
                  x1[nu] += 1; // x - mu + nu
                  x2[mu] -= 1; // x - mu
                  leaf += U(x, nu) * U(x1, mu).dagger() * U(x2, nu).dagger() * U(x2, mu);
                  x1[mu] += 1; // x + nu
                  x1[nu] -= 1; // x
                  x2[mu] += 1; // x

                  x1[mu] -= 1; // x - mu
                  x2[mu] -= 1; // x - mu
                  x2[nu] -= 1; // x - mu -nu
                  x3[nu] -= 1; // x - nu
                  leaf += U(x1, mu).dagger() * U(x2, nu).dagger() * U(x2, mu) * U(x3, nu);
                  x1[mu] += 1; // x
                  x2[mu] += 1; // x + nu
                  x2[nu] += 1; // x
                  x3[nu] += 1; // x

                  x1[nu] -= 1; // x - nu
                  x2[nu] -= 1; // x - nu
                  x2[mu] += 1; // x + mu -nu
                  leaf += U(x1, nu).dagger() * U(x1, mu) * U(x2, nu) * U(x, mu).dagger();
                  x1[nu] += 1; // x
                  x2[nu] += 1; // x + mu
                  x2[mu] -= 1; // x
                }
                // traceless and anti-hermitian
                // here we include a factor 1/2 already
                G[mu][nu] = traceless_antiherm(leaf);
                // trace(G_{mu,nu}^a G_{mu,nu}^a)
                res += retrace(G[mu][nu] * G[mu][nu]);
              }
            }

            if (U.getndims() == 4) {
              // sum up the topological charge contribution now
              for (int i = 0; i < eps4.N; i++) {
                int i1 = eps4.eps_idx[i][0];
                int i2 = eps4.eps_idx[i][1];
                int i3 = eps4.eps_idx[i][2];
                int i4 = eps4.eps_idx[i][3];

                // when Gmunu components from the lower triangle are to be used,
                // we can simply skip them and multiply our normalisation by a factor of
                // four in total
                if (i2 < i1) {
                  continue;
                }
                if (i4 < i3) {
                  continue;
                }
                Q += eps4.eps_val[i] * retrace(G[i1][i2] * G[i3][i4]);
              }
            }
            if (U.getndims() == 2) {
              Q += -std::imag(trace((G[0][1] - G[1][0])));
            }
          }
        }
      }
    }
    // now we need to devide by 2, but we get a factor of two since we only
    // averaged mu < nu
    res = -res / U.getVolume();
    // averaged over four plaquette Wilson loops 1./4./4.
    if (cloverdef) {
      res /= 16.;
    }
    // factor 4 from summing only mu < nu and rho < sigma
    Q = -4. * Q / (32.0 * M_PI * M_PI);
    // factor 1/16 from G_\mu\nu with clover definition
    if (cloverdef)
      Q /= 16.;
  }

} // namespace flat_spacetime
