/**
 * @file parameters.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief global parameters of the programs
 * This file contains a set of 'struct' which serve as containers for the parameters of
 * the programs.
 * @version 0.1
 * @date 2022-05-30
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <string>

namespace global_parameters {

  // default string width for the beta name in output file
  const size_t g_beta_str_width = 6;

  /* struct for parameters concerning the physics of the system */
  struct physics {
  public:
    size_t Lx; // spatial  lattice size X > 0
    size_t Ly; // spatial  lattice size Y > 0
    size_t Lz; // spatial  lattice size Z > 0
    size_t Lt; // temporal lattice size T > 0
    size_t ndims = 4; // number of dimensions, 2 <= ndims <= 4

    bool flat_metric = true; // false when considering spacetime curvature
    bool rotating_frame = false; // true when we consider a rotating lattice
    double Omega = 0.0; // imaginary angular frequency (rotating lattice)

    bool include_gauge = false; // true when the action contains the gauge part
    double beta; // beta value

    bool include_staggered_fermions =
      false; // true when the action contains staggered fermions.
    double m0; // bare quark mass
    double xi = 1.0; // anisotropy parameter
    bool anisotropic = false; // use anisotropic lattice
  };

  struct measure_glueball_u1 {
    bool do_measure = false; // measure the glueball correlator
    bool doAPEsmear = false; // apply APE smearing to the links
    size_t nAPEsmear; // number of APE smearing iterations
    double alphaAPEsmear; // alpha parameter for the smearing. alpha=1 -> no smearing.
    size_t max_length_loops = 4; // maximum length of loops interpolating the glueballs
    std::string subfld = "./"; // subfolder for the measure of glueball correlators only
    bool lenghty_conf_name = true; // add measure information in correlator filename
  };

  /* optional parameters for the measure program the in U(1) theory */
  struct measure_u1 {
    size_t verbosity = 0; // verbosity of the output

    size_t n_meas = 10; // total number of sweeps
    size_t icounter = 0; // initial counter for updates
    size_t seed = 13526463; // PRNG seed

    size_t nstep = 1; // measure each nstep config
    bool Wloop = false; // whether to measure the Wilson loops or not
    bool gradient_flow = false; // whether to measure the gradient flow or not
    double epsilon_gradient_flow =
      0.01; // step size in the integration of the gradient flow equations
    double tmax = 1.0; // tmax for gradient flow
    std::string conf_dir = "./"; // directory where gauge configurations are stored
    std::string res_dir = "./"; // directory where results from measurements for
                                // potential, potentialsmall are stored

    bool potential =
      false; // measure potential: the loops W(x, t, y=z=0) and W(x, y, t=z=0) are
             // measured with a maximum size of lattice extent * sizeloops, and written to
             // separate files. Only available for ndims=3,4
    bool potentialsmall =
      false; // The loops W(x, t, y) are measured up to x, y=min(4, lattice extent), t <=
             // Lt * sizeloops and saved to one file. Only available for ndim=3
    bool append = false; // are measurements for potential appended to an existing file,
                         // or should it be overwritten?
    double sizeWloops =
      0.5; // Wilson-Loops are measured up to this fraction of the lattice extent
    size_t n_apesmear =
      0; // number of APE smearings done on the lattice before measurement.
    // APE-smearing is done before measuring the potential and small potential, it does
    // not affect the gradient flow and Wilson-loops
    double alpha = 1.0; // parameter alpha for APE smearings
    bool smear_spatial_only = false; // should smearing be done only for spacial links?
    bool smear_temporal_only = false; // should smearing be done only for temporal links?

    std::string conf_basename = "u1_conf"; // root of the output files names
    bool lenghty_conf_name = true; // add ensemble information in configuration name
    size_t beta_str_width = g_beta_str_width; // length of the beta value config filename

    bool pion_staggered =
      false; // whether to measure the staggered pion correlator or not
    double m0; // bare quark mass
    std::string solver = "CG"; // Type of solver: CG, BiCGStab
    double tolerance_cg = 1e-10; // Tolerance for the solver for the dirac operator
    size_t solver_verbosity = 0; // Verbosity for the solver for the dirac operator
    size_t seed_pf = 97234719; // Seed for the evaluation of the fermion determinant

    measure_glueball_u1
      measure_glueball_params; // structure for the measure of the glueball
  };

  /* Optional parameters for the hmc the in U(1) theory */
  struct hmc_u1 {
    bool do_hmc = true; // whether to do the hmc evolution or not

    size_t N_save = 100; // N_save
    size_t n_meas = 10; // total number of sweeps
    size_t icounter = 0; // initial counter for updates
    size_t seed = 13526463; // PRNG seed
    bool heat = true; // randomness of the initial config, true: hot, false: cold
    bool restart =
      false; // restart from the last configuration (reads from nconf_counter.txt)
    bool acceptreject = true; // use accept-reject
    std::string configfilename = ""; // configuration filename used in case of restart

    size_t N_rev = 0; // frequency of reversibility tests N_rev, 0: no reversibility test
    size_t n_steps = 1000; // n_steps
    double tau = 1; // trajectory length tau
    size_t exponent = 0; // exponent for rounding
    size_t integs = 0; // itegration scheme to be used: 0=leapfrog, 1=lp_leapfrog, 2=omf4,
                       // 3=lp_omf4, 4=Euler, 5=RUTH, 6=omf2

    // itegration scheme to be used:
    // leapfrog, lp_leapfrog, omf4, lp_omf4, Euler, RUTH, omf2
    std::string integrator = "leapfrog";

    std::string solver = "CG"; // Type of solver: CG, BiCGStab
    double tolerance_cg = 1e-10; // Tolerance for the solver for the dirac operator
    size_t solver_verbosity = 0; // Verbosity for the solver for the dirac operator
    size_t seed_pf = 97234719; // Seed for the evaluation of the fermion determinant
    std::string conf_dir = "."; // Output directory

    std::string conf_basename = "u1_conf"; // root of the output files names
    bool lenghty_conf_name = true; // add ensemble information in configuration name
    size_t beta_str_width = g_beta_str_width; // length of the beta value config filename

    bool make_omeas = false; // true/false when online measurement are ON/OFF
    measure_u1 omeas; // stuct for online measurements
  };

  /* optional parameters for the MCMC the in U(1) theory */
  struct metropolis_u1 {
    size_t n_meas = 10; // total number of sweeps
    size_t icounter = 0; // initial counter for updates
    size_t seed = 13526463; // PRNG seed

    size_t N_save = 100; // save each N_save config

    bool restart = false; // restart from an existing configuration
    std::string configfilename = ""; // configuration filename used in case of restart
    std::string conf_dir = "./"; // directory where gauge configurations are stored
    std::string conf_basename = "u1_conf"; // root of the output files names
    bool lenghty_conf_name = true; // add ensemble information in configuration name
    size_t beta_str_width = g_beta_str_width; // length of the beta value config filename

    size_t N_hit = 10; // N_hit updates are performed on each link during one sweep
    double heat =
      1.0; // determines if thermalization starts from a hot (=1) or cold(=0) config
    double delta =
      1.0; // quantifies how much the prooposed new link can differ from the current link
  };

} // namespace global_parameters
