// parameters.hh
/**
 * @brief global parameters of the programs
 * This file contains a set of 'struct' which serve as containers for the global
 * parameters of the programs.
 */

#pragma once

#include <string>

namespace global_parameters {

  /* struct for general parameters common to all programs */
  struct general {
  public:
    size_t Lx; // spatial lattice size x > 0
    size_t Ly; // spatial lattice size y > 0
    size_t Lz; // spatial lattice size z > 0
    size_t Lt; // temporal lattice size > 0
    size_t ndims; // number of dimensions, 2 <= ndims <= 4
    size_t N_save = 1000; // N_save
    size_t N_meas = 10; // total number of sweeps
    size_t icounter = 0; // initial counter for updates
    double beta; // beta value
    double m0; // bare quark mass
    size_t seed = 13526463; // PRNG seed
    double heat = 1.0; // randomness of the initial config, 1: hot, 0: cold
    bool restart = false; // restart from an existing configuration
    bool acceptreject = true; // use accept-reject
    std::string configfilename; // configuration filename used in case of restart
  };

  /* optional parameters for the hmc the in U(1) theory */
  struct hmc_u1 {

    size_t N_rev = 0; // frequency of reversibility tests N_rev, 0: no reversibility test
    size_t n_steps = 1000; // n_steps
    double tau = 1; // trajectory length tau
    size_t exponent = 0; //exponent for rounding
    size_t integs = 0; // itegration scheme to be used: 0=leapfrog, 1=lp_leapfrog, 2=omf4, 3=lp_omf4, 4=Euler, 5=RUTH, 6=omf2
    bool no_fermions = 0; // Bool flag indicating if we're ignoring the fermionic action.
    std::string solver = "CG"; // Type of solver: CG, BiCGStab
    double tolerance_cg = 1e-10; // Tolerance for the solver for the dirac operator
    size_t solver_verbosity = 0; // Verbosity for the solver for the dirac operator
    size_t seed_pf = 97234719; // Seed for the evaluation of the fermion determinant
    std::string outdir = "."; // Output directory
  };

} // namespace global_parameters