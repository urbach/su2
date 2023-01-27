#pragma once
// #include"gradient_flow.hh"
#include "adjointfield.hh"
#include "flat-energy_density.hh"
#include "gauge_energy.hpp"
#include "gaugeconfig.hh"
#include "gaugemonomial.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "su2.hh"
#include "update_gauge.hh"

#include <fstream>
#include <iomanip>
#include <string>

template <typename Float, class Group>
void runge_kutta(hamiltonian_field<Float, Group> &h,
                 monomial<Float, Group> &SW,
                 const double eps) {
  double zfac[5] = {(-17.0) / (36.0), (8.0) / (9.0), (-3.0) / (4.0)};
  double expfac[3] = {-36.0 / 4. / 17.0, 1., -1.};

  // following arxiv:1006.4518
  //
  // W0      = V_t
  // W1      = exp(1/4 Z0) W0
  // W2      = exp(8/9 Z1 - 17/36 Z0) W1 = exp(Z') W1
  // V_t+eps = exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2 = exp(3/4 Z2 - Z') W2
  //
  // Zi = eps*Z(Wi)
  // before the three steps zero the derivative field
  zeroadjointfield(*(h.momenta));

  const double fact_f = 2.0 * h.Fact_Nc_force_Luscher * h.U->getNc() / h.U->getBeta();
  for (int f = 0; f < 3; f++) {
    // evolution of eq. C.2 of https://arxiv.org/pdf/1006.4518.pdf
    // we have to cancel beta/N_c from the derivative

    // updating the momenta *(h.momenta)
    SW.derivative(*(h.momenta), h, fact_f * zfac[f]);
    // update the flowed gauge field Vt
    update_gauge(h, -eps * expfac[f]);
  }
  return;
}

namespace flat_spacetime {

  /**
   * @brief Printing Wilson gradient flow evolution on 'path'
   * Prints a dataframe for the Wilson flow evolution of
   * flowtime, Plaquette(full, spatial-spatial, temporaVt-save(path+".conf"); ial),
   * Energy_plaquette(...), Energy_clover[improved formula](...), etc. Notes:
   * - Each type of plaquette has limit 1 --> different normalization factors for P_ss and
   * P_ts.
   * - E_i = 1 - P_i, i="","ss","ts"
   * - E = E_ss + E_ts
   * @tparam Group
   * @param U
   * @param path
   * @param tstart 1st value of flow time (when loading flowed gauge configuration)
   * @param tmax
   * @param eps
   */
  template <class Group>
  void gradient_flow(const gaugeconfig<Group> &U,
                     std::string const &path,
                     const double &tmax,
                     const double &eps,
                     const double &xi,
                     const double &tstart,
                     const bool &save_config) {
    const size_t d = U.getndims();
    const double ndims_fact = spacetime_lattice::num_pLloops_half(d);
    const double ndims_fact_ss = spacetime_lattice::num_pLloops_half(d - 1);

    double t[3];
    double P[3], E[3], Q[3];
    double P_ss[3], E_ss[3], Q_ss[3]; // spatial-spatial contribution

    std::ostringstream oss;

    double density = 0., topQ = 0.; // dummy variables
    double density_ss = 0., topQ_ss = 0.; // dummy variables
    for (unsigned int i = 0; i < 3; i++) {
      t[i] = tstart;
      P[i] = 0.;
      E[i] = 0.;
      Q[i] = 0.;
      P_ss[i] = 0.;
      E_ss[i] = 0.;
      Q_ss[i] = 0.;
    }

    const double den_common = U.getVolume() * double(U.getNc());

    P[2] = flat_spacetime::gauge_energy(U) / den_common;
    P_ss[2] = flat_spacetime::gauge_energy(U, true) / den_common;
    flat_spacetime::energy_density(U, density, topQ);
    E[2] = density;
    flat_spacetime::energy_density(U, density_ss, topQ_ss, true, true);
    E_ss[2] = density_ss;

    // definine a fictitious gauge configuration Vt, momenta and hamiltonian field to
    // evolve with the flow
    gaugeconfig<Group> Vt(U);
    adjointfield<double, Group> deriv(U.getLx(), U.getLy(), U.getLz(), U.getLt(), d);
    hamiltonian_field<double, Group> h(deriv, Vt);
    gaugemonomial<double, Group> SW(0, xi); // Wilson (pure) gauge action

    if (tstart == 0.0) {
      oss << "t "; // flow time
      oss << "xi "; // e_ts(t)/e_ss(t)
      oss << "P P_ss "; // plaquette
      oss << "Ep Ep_ss "; // energy from regular plaquette
      oss << "Ec Ec_ss "; // energy from clover-leaf plaquette
      oss << "Q" << std::endl;
    }

    // evolution of t[1] until tmax
    //(note: eps=0.01 and tmax>0 --> the loop ends at some point)
    // at each step we consider a triplet of values for t,P,E,Q

    while (t[1] < tmax - eps) {
      // splicing the results at t[2] to the new 0-th temporal time slice
      t[0] = t[2];
      P[0] = P[2];
      E[0] = E[2];
      Q[0] = Q[2];
      P_ss[0] = P_ss[2];
      E_ss[0] = E_ss[2];

      for (unsigned int x0 = 1; x0 < 3; x0++) {
        t[x0] = t[x0 - 1] + eps;
        runge_kutta(h, SW, eps); // apply Runge-Kutta integration method
        P[x0] = flat_spacetime::gauge_energy(Vt) / den_common;
        P_ss[x0] = flat_spacetime::gauge_energy(Vt, true) / den_common;
        flat_spacetime::energy_density(Vt, density, topQ, true);
        E[x0] = density;
        Q[x0] = topQ;
        flat_spacetime::energy_density(Vt, density_ss, topQ_ss, true, true);
        E_ss[x0] = density_ss;
      }

      const double tsqr = t[1] * t[1];

      const double Ep = ndims_fact - P[1];
      const double Ep_ss = ndims_fact_ss - P_ss[1];

      const double Ep_ts = Ep - Ep_ss;

      const double Ep_ts_bar = Ep_ts / (d - 1);
      const double Ep_ss_bar = Ep_ss / ndims_fact_ss;

      const double xi2_R = Ep_ts_bar / Ep_ss_bar;
      const double xi_R = sqrt(xi2_R);

      const double t2Ep = tsqr * Ep;
      const double t2Ep_ss = tsqr * Ep_ss;
      const double t2Ep_ts = tsqr * Ep_ts;

      const double Ec = E[1];
      const double Ec_ss = E_ss[1];
      const double Ec_ts = Ec - Ec_ss;

      const double t2Ec = tsqr * Ec;
      const double t2Ec_ss = tsqr * Ec_ss;

      oss << std::scientific; // using scientific notation
      oss.precision(16);
      oss << t[1] << " ";
      oss << xi_R << " ";
      oss << P[1] << " " << P_ss[1] << " ";
      oss << Ep << " " << Ep_ss << " ";
      oss << Ec << " " << Ec_ss << " ";
      oss << Q[1] << "\n";
    }

    std::ofstream ofs;
    if (tstart == 0.0) {
      ofs.open(path, std::ios::out);
    } else {
      ofs.open(path, std::ios::app);
    }

    ofs << oss.str();
    if (save_config) {
      Vt.save(path + "_t" + std::to_string(tmax) + ".conf");
    }

    return;
  }

} // namespace flat_spacetime
