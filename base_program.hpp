/**
 * @file base_program.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief mother class for all the classes that handle the programs runs
 * @version 0.1
 * @date 2022-09-01
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include "flat-energy_density.hh"
#include "gauge_energy.hh"
#include "flat-sweep.hh" // flat spacetime
#include "gaugeconfig.hh"
#include "io.hh"
#include "omeasurements.hpp"
#include "parse_input_file.hh"
#include "random_gauge_trafo.hh"
// #include "rotating-energy_density.hpp" // rotating spacetime
// #include "rotating-gauge_energy.hpp" // rotating spacetime
// #include "rotating-sweep.hpp" // rotating spacetime
#include "su2.hh"
#include "u1.hh"
#include "vectorfunctions.hh"
#include "version.hh"

#ifdef _USE_OMP_
#include <omp.h>
#endif

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>

namespace po = boost::program_options;

/**
 * @brief parsing the command line for the main program.
 * Parameters available: "help" and "file"
 * see implementation for details
 * @param ac argc from the standard main() function
 * @param av argv from the standard main() function
 */
void parse_command_line(int ac, char *av[], std::string &input_file) {
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce this help message")(
    "file,f", po::value<std::string>(&input_file)->default_value("NONE"),
    "yaml input file");

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(0);
  }
  return;
}

/**
 * @brief struct of boolean flags
 * each flag is true when the corresponding program is running.
 * NOTE: Only one flat at the time can be true
 */
struct running_program {
  bool do_hmc = false; // Hybrid Monte Carlo algorithm
  bool do_metropolis = false; // Metropolis algorithm
  bool do_omeas = false; // only measuring observables
};

/**
 * @brief node without conflicting nodes
 * In the input file one may have specified multiple MC algorithms or none (only offline
 * measurements). This function takes the path to theinput file as input and returns a
 * YAML node cleaned and ready to be used.
 *
 * @param rp running program
 * @param input_file path to the input file
 */
YAML::Node get_cleaned_input_file(running_program &rp, const std::string &input_file) {
  YAML::Node nd = YAML::LoadFile(input_file);
  YAML_parsing::inspect_node in(nd);

  bool &do_hmc = rp.do_hmc;
  bool &do_metropolis = rp.do_metropolis;
  bool &do_omeas = rp.do_omeas;

  if (nd["hmc"]) {
    in.read_verb<bool>(do_hmc, {"hmc", "do_mcmc"});
    if (!do_hmc) {
      nd.remove("integrator");
      nd.remove("hmc");
    }
  }

  do_metropolis = false;
  if (nd["metropolis"]) {
    in.read_verb<bool>(do_metropolis, {"metropolis", "do_mcmc"});
    if (!do_metropolis) {
      nd.remove("metropolis");
    }
  }

  do_omeas = bool(nd["omeas"]);
  std::string err = ""; // error message
  try {
    if (do_hmc && do_metropolis) { // both options are incompatible
      err = "ERROR: Can't run simultaneously hmc and metropolis algorithms.\n";
      throw err;
    } else if (!do_hmc && !do_metropolis && !nd["omeas"]["offline"]) {
      err = "Error: You must write the 'offline' measurements node inside 'omeas' "
            "because you're not running any MCMC algorithm.\n";
      throw err;
    }
  } catch (...) {
    std::cerr << err << "Check your input file: " << input_file << "\n";
    std::cerr << "Aborting.\n";
    std::abort();
  }

  return nd;
}

namespace gp = global_parameters;

template <class Group, class sparam_type> class base_program {
protected:
  std::string algo_name =
    "UNNAMED_PROGRAM"; // name of the algorithm; initialized specified in child classes

  gp::physics pparams; // physics parameters
  sparam_type sparams; // specific parameters to the given run
  gp::measure omeas; // omeasurements parameters
  //    YAML::Node nd; // yaml node

  size_t threads;

  // filename needed for saving results from potential and potentialsmall
  std::string filename_fine;
  std::string filename_coarse;
  std::string filename_nonplanar;

  std::string conf_path_basename; // basename for configurations
  gaugeconfig<Group> U; // gauge configuration evolved through the algorithm

  bool g_heat; // hot or cold starting configuration
  size_t g_icounter = 0; // 1st configuration(trajectory) to load from

  double normalisation;
  size_t facnorm;

  std::ofstream os;
  std::ofstream acceptancerates;

public:
  base_program() {}
  ~base_program() {}

  virtual void print_program_info() const = 0;

  void print_git_info() const {
    std::cout << "## GIT branch " << GIT_BRANCH;
    std::cout << " on commit " << GIT_COMMIT_HASH << "\n";
  }

  void print_info() const {
    this->print_program_info();
    this->print_git_info();
  }

  virtual void parse_input_file(const YAML::Node &nd) = 0;

  /**
   * @brief create the necessary output directories
   */
  void create_directories() {
    namespace fsys = boost::filesystem;
    fsys::create_directories(fsys::absolute(sparams.conf_dir));
    fsys::create_directories(fsys::absolute(omeas.res_dir));
  }

  void set_omp_threads() {
#ifdef _USE_OMP_
    /**
     * the parallelisation of the sweep-function first iterates over all odd points in t
     * and then over all even points because the nearest neighbours must not change
     * during the updates, this is not possible for an uneven number of points in T
     * */
    if (pparams.Lt % 2 != 0) {
      std::cerr << "For parallel computing an even number of points in T is needed!"
                << std::endl;
      omp_set_num_threads(1);
      std::cerr << "Continuing with one thread." << std::endl;
    }
    // set things up for parallel computing in sweep
    threads = omp_get_max_threads();
#else
    threads = 1;
#endif
    std::cout << "threads " << threads << std::endl;
  }

  /**
   * @brief initialize the potential filenames attributes
   */
  void set_potential_filenames() {
    // filename needed for saving results from potential and potentialsmall
    filename_fine = io::measure::get_filename_fine(pparams, omeas);
    filename_coarse = io::measure::get_filename_coarse(pparams, omeas);
    filename_nonplanar = io::measure::get_filename_nonplanar(pparams, omeas);

    // write explanatory headers into result-files, also check if measuring routine is
    // implemented for given dimension
    if (omeas.potentialplanar) {
      io::measure::set_header_planar(pparams, omeas, filename_coarse, filename_fine);
    }
    if (omeas.potentialnonplanar) {
      io::measure::set_header_nonplanar(pparams, omeas, filename_nonplanar);
    }
    return;
  }

  double gauge_energy(const gp::physics &pparams,
                      const gaugeconfig<Group> &U,
                      const bool spatial_only = false) {
    if (pparams.flat_metric) {
      return flat_spacetime::gauge_energy(U, spatial_only);
    }
    if (pparams.rotating_frame) {
      spacetime_lattice::fatal_error("Rotating metric not supported yet.", __func__);
      // return rotating_spacetime::gauge_energy(U, pparams.Omega, spatial_only);
    } else {
      spacetime_lattice::fatal_error("Invalid metric when calling: ", __func__);
      return {};
    }
  }

  template <class T>
  void energy_density(const gp::physics &pparams,
                      const gaugeconfig<T> &U,
                      double &E,
                      double &Q,
                      const bool &cloverdef = true,
                      const bool &ss = false) {
    if (pparams.flat_metric) {
      flat_spacetime::energy_density(U, E, Q, cloverdef, ss);
    }
    if (pparams.rotating_frame) {
      spacetime_lattice::fatal_error("Rotating metric not supported yet.", __func__);
      // rotating_spacetime::energy_density(U, pparams.Omega, E, Q, cloverdef);
    }
    return;
  }

  /**
   * @brief create gauge configuration with correct geometry
   */
  void create_gauge_conf() {
    gaugeconfig<Group> U0(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims,
                          pparams.beta);
    U = U0;
  }

  /**
   * @brief initialize the gauge configuration for the Markov chain Monte Carlo programs
   */
  void init_gauge_conf_mcmc() {
    /**
     * @brief measuring spatial plaquettes only means only (ndims-1)/ndims of all
     * plaquettes are measured, so need facnorm for normalization to 1
     *
     */
    facnorm = (pparams.ndims > 2) ? pparams.ndims / (pparams.ndims - 2) : 0;

    if (sparams.restart) {
      std::cout << "## restart " << sparams.restart << "\n";
      const std::vector<std::string> v_ncc = io::read_nconf_counter(sparams.conf_dir);
      g_heat = boost::lexical_cast<bool>(v_ncc[0]);
      g_icounter = std::stoi(v_ncc[1]);
      std::string config_path = v_ncc[2];

      const size_t err = U.load(config_path);
      if (err != 0 && sparams.do_mcmc) {
        std::cout << "Error: failed to load initial gauge configuration for "
                     "intializing the Markov chain Monte Carlo. Aborting.\n";
        std::abort();
      }
    } else {
      g_heat = (sparams.heat == true) ? 1.0 : 0.0;
      g_icounter = 0;
      hotstart(U, sparams.seed, g_heat);
    }

    double plaquette = flat_spacetime::gauge_energy(U);
    double fac = 2. / U.getndims() / (U.getndims() - 1);
    normalisation = fac / U.getVolume() / double(U.getNc());

    std::cout << "## Normalization factor: A = 2/(d*(d-1)*N_lat*N_c) = "
              << std::scientific << std::setw(18) << std::setprecision(15)
              << normalisation << "\n";
    std::cout << "## Acceptance rate parcentage: rho = rate/(i+1)\n";

    std::cout << "## Initial Plaquette: " << plaquette * normalisation << std::endl;

    random_gauge_trafo(U, 654321);
    plaquette = flat_spacetime::gauge_energy(U);
    std::cout << "## Plaquette after rnd trafo: " << plaquette * normalisation
              << std::endl;
  }

  void open_output_data() {
    // doing only offline measurements
    if (!sparams.do_mcmc) {
      return;
    }

    const std::string file = sparams.conf_dir + "/output.u1-" + algo_name + ".data";
    if (g_icounter == 0) {
      os.open(file, std::ios::out);
    } else {
      os.open(file, std::ios::app);
    }

    return;
  }

  /**
   * @brief part of the program flow common to all programs
   */
  void pre_run(const YAML::Node &nd) {
    this->print_program_info();
    this->print_git_info();

    // printing the yaml main node -> reproducibility of the run
    std::cout << "## Cleaned yaml node:\n";
    std::cout << nd << "\n";

    this->parse_input_file(nd);
    this->create_gauge_conf();

    this->create_directories();

    this->open_output_data();
  }

  /**
   * @brief run the (generic) program
   *
   * @param path path to the input file
   */
  virtual void run(const YAML::Node &nd) = 0;

  /**
   * @brief online measurements over the i-th trajectory
   *
   * @param i trajectory index
   */
  void do_omeas_i(const size_t &i) {
    const bool flag_i =
      sparams.do_omeas && (i > omeas.icounter) && ((i % omeas.nstep) == 0);

    if (!flag_i) {
      return;
    }

    if (!sparams.do_mcmc) { // doing only offline measurements
      const std::string path_i = conf_path_basename + "." + std::to_string(i);
      int ierrU = U.load(path_i);

      if (ierrU == 1) { // cannot load gauge config
        return; // simply ignore configuration
      }
    }

    if (omeas.potentialplanar || omeas.potentialnonplanar) {
      gaugeconfig<Group> U1 = U;
      // smear lattice
      for (size_t smears = 0; smears < omeas.n_apesmear; smears += 1) {
        APEsmearing<double, Group>(U1, omeas.alpha, omeas.smear_spatial_only);
      }
      if (omeas.potentialplanar) {
        omeasurements::meas_loops_planar_pot(U1, pparams, omeas.sizeWloops,
                                             (*this).filename_coarse,
                                             (*this).filename_fine, i);
      }

      if (omeas.potentialnonplanar) {
        omeasurements::meas_loops_nonplanar_pot(U1, pparams, omeas.sizeWloops,
                                                (*this).filename_nonplanar, i);
      }
    }

    if ((*this).omeas.Wloop) {
      if ((*this).omeas.verbosity > 0) {
        std::cout << "## online measuring: Wilson loop\n";
      }
      omeasurements::meas_wilson_loop<Group>(U, i, omeas.res_dir);
    }
    if ((*this).omeas.gradient_flow.measure_it) {
      if ((*this).omeas.verbosity > 0) {
        std::cout << "## online measuring: Gradient flow\n";
      }
      omeasurements::meas_gradient_flow<Group>(U, i, pparams, (*this).omeas);
    }

    if ((*this).omeas.pion_staggered) {
      if ((*this).omeas.verbosity > 0) {
        std::cout << "## online measuring: Pion correlator\n";
      }
      omeasurements::meas_pion_correlator<Group>(U, i, pparams.m0, (*this).omeas);
    }

    if ((*this).omeas.glueball.do_measure) {
      if ((*this).omeas.verbosity > 0) {
        std::cout << "## online measuring: J^{PC} glueball correlators.\n";
      }

      const std::string glb_interp = omeas.glueball.interpolator_type;
      omeasurements::meas_glueball_correlator<Group>(glb_interp, U, i, (*this).omeas);
    }
    return;
  }

  /**
   * @brief operations to be done after the i-th step of the MCMC
   *
   * @param i configuration index
   */
  void after_MCMC_step(const size_t &i, const bool &do_omeas) {
    if (i > 0 && (i % (*this).sparams.N_save) ==
                   0) { // saving (*this).U after each N_save trajectories
      std::string path_i = (*this).conf_path_basename + "." + std::to_string(i);
      if ((*this).sparams.do_mcmc) {
        (*this).U.save(path_i);
      }

      if (do_omeas) {
        this->do_omeas_i(i);
      }

      if ((*this).sparams.do_mcmc) {
        // storing last conf index (only after online measurements has been done)
        io::update_nconf_counter((*this).sparams.conf_dir, (*this).g_heat, i, path_i);
      }
    }
  }
};
