// parse_input_file.cc
/**
 * @file parse_input_file.cc
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief definitions of parse_input_file.hh
 * @version 0.1
 * @date 2022-05-04
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <iostream>
#include <string>

#include "geometry.hh"
#include "parse_input_file.hh"

namespace input_file_parsing {

  void parse_command_line(int ac, char *av[], std::string &input_file) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "produce this help message")(
      "file,f", po::value<std::string>(&input_file)->default_value("NONE"),
      "yaml input file");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      std::abort();
    }

    if (input_file == "NONE") {
      std::cout << "Error: invalid command line parameters.\nCheck using the "
                   "\"-h/--help\" option."
                << "\n";
      std::abort();
    }

    return;
  }

  int validate_geometry(gp::physics &pparams) {
    if (pparams.ndims > 4 || pparams.ndims < 2) {
      std::cerr << "2 <= ndims <= 4!" << std::endl;
      return 1;
    }
    if (pparams.Lx < 1 || pparams.Ly < 1 || pparams.Lz < 1 || pparams.Lt < 1) {
      std::cerr << "All box extents must be > 1!" << std::endl;
      return 1;
    }
    if (pparams.ndims < 4) {
      std::cerr << "## Warning: ndims==" << pparams.ndims << " --> flattening the \'z\' ";
      pparams.Lz = 1;
      std::string s_end = "";
      if (pparams.ndims < 3) {
        pparams.Ly = 1;
        std::cerr << "and \'y\' ";
        s_end = "s";
      }
      std::cerr << "direction" << s_end << "\n";
    }

    return 0;
  }

  /**
   * @brief check if beta_str_width >= 1, otherwise it aborts
   * @param n width (number of digits after the dot)
   */
  void validate_beta_str_width(const size_t &n) {
    if (n < 1) {
      std::cerr << "Error: beta string width in the output file should have at least "
                   "decimal 1 digit.";
      std::cerr << "Aborting.\n";
      std::abort();
    }
    return;
  }

  /**
   * @brief check if N_hit >= 1, otherwise it aborts
   * @param n N_hit (number of times an update per link is attempted)
   */
  void validate_N_hit(const size_t &n) {
    if (n < 1) {
      std::cerr << "Error: N_hit should be at least 1, otherwise nothing will happen";
      std::cerr << "Aborting.\n";
      std::abort();
    }
    return;
  }

  namespace Yp = YAML_parsing;

  void parse_geometry(Yp::inspect_node &in, gp::physics &pparams) {
    YAML::Node R = in.get_root();

    const bool spec_L = (bool)R["geometry"]["L"];
    const bool spec_Lxyz =
      (R["geometry"]["X"] || R["geometry"]["Y"] || R["geometry"]["Z"]);

    if (spec_L ^ spec_Lxyz) {
      // either L^3*T or Lx*Ly*Lz*T have been specified (but not both)

      if (spec_L) {
        size_t L;
        in.read_verb<size_t>(L, {"geometry", "L"});
        pparams.Lx = L;
        pparams.Ly = L;
        pparams.Lz = L;
      } else {
        in.read_verb<size_t>(pparams.Lx, {"geometry", "X"});
        in.read_verb<size_t>(pparams.Ly, {"geometry", "Y"});
        in.read_verb<size_t>(pparams.Lz, {"geometry", "Z"});
      }

    } else {
      std::cerr << "Error: check your input file. ";
      std::cerr
        << "Either you specify L=Lx=Ly=Lz or each dimension separately, not both.\n";
      std::cerr << "Aborting.\n";
      std::abort();
    }

    in.read_verb<size_t>(pparams.Lt, {"geometry", "T"});
    in.read_verb<size_t>(pparams.ndims, {"geometry", "ndims"});

    int gerr = validate_geometry(pparams);
    if (gerr > 0) {
      std::cerr << "Error: invalid geometry parameters. Check L (or X,Y,Z), T, ndims "
                   "in your input file.";
      std::cerr << "Aborting.\n";
      std::abort();
    }

    if (R["geometry"]["rotating_frame"]) {
      pparams.rotating_frame = true;
      pparams.flat_metric = false; // metric is not flat
      in.read_verb<double>(pparams.Omega, {"geometry", "rotating_frame", "Omega"});
    }

    return;
  }

  /**
   * @brief checks that `heat` and `restart` are not both present at the same time
   *
   */
  void heat_restart_incompatible(const YAML::Node &nd) {
    if (nd["restart"]) {
      if (!(nd["restart"].as<bool>()) && !nd["heat"]) {
        std::cerr << "Error: Your gave 'restart==false' but didn't specify the 'heat' "
                     "parameter.\n";
        std::cerr << "Aborting.\n";
        std::abort();
      }
      if (nd["restart"] && nd["heat"]) {
        std::cerr << "Error: "
                  << "'restart' and 'heat' conditions are incompatible in the hmc. "
                  << "Aborting.\n";
        std::abort();
      }
    }
    if (!nd["restart"] && !nd["heat"]) {
      std::cerr << "Error: "
                << "Please pass either 'restart' or 'heat' to the hmc. "
                << "Aborting.\n";
      std::abort();
    }
  }

  /**
   * @brief parsing parameters from node of glueball correlator
   * @param in input file parser
   * @param mgparams measure_glueball_parameters structure
   */
  void parse_glueball_measure(Yp::inspect_node &in,
                              const std::vector<std::string> &inner_tree,
                              gp::measure_glueball &mgparams) {
    const std::vector<std::string> state0 = in.get_InnerTree();
    in.dig_deeper(inner_tree); // entering the glueball node
    YAML::Node nd = in.get_outer_node();

    mgparams.interpolator = bool(nd["interpolator"]);
    if (mgparams.interpolator) {
      in.read_verb<std::string>(mgparams.interpolator_type, {"interpolator", "type"});
      in.read_verb<bool>(mgparams.spatial_loops, {"interpolator", "spatial"});
      in.read_verb<size_t>(mgparams.rmin, {"interpolator", "rmin"});
      in.read_verb<size_t>(mgparams.rmax, {"interpolator", "rmax"});
      in.read_verb<bool>(mgparams.save_interpolator, {"interpolator", "save"});
      in.read_verb<bool>(mgparams.correlator, {"correlator"});
    } else {
      std::cerr << "Error: No 'glueball:interpolator' node found in the input file.\n";
      std::cerr << "Aborting.\n";
      std::abort();
    }

    in.read_verb<bool>(mgparams.correlator, {"correlator"});
    mgparams.do_measure =
      mgparams.interpolator || mgparams.correlator; // true if we save one of the two
    in.read_verb<bool>(mgparams.doAPEsmear, {"do_APE_smearing"});
    if (mgparams.doAPEsmear) {
      in.read_sequence_verb<size_t>(mgparams.vec_nAPEsmear, {"APE_smearing", "n"});
      in.read_verb<double>(mgparams.alphaAPEsmear, {"APE_smearing", "alpha"});
    }
    in.read_opt_verb<bool>(mgparams.lengthy_file_name, {"lengthy_file_name"});
    in.read_opt_verb<bool>(mgparams.use_res_dir, {"res_dir"});

    in.set_InnerTree(state0); // reset to previous state
  }

  /**
   * @brief parsing parameters from node of gradient flow
   * @param in input file parser
   * @param mgparams measure_gradient_flow structure
   */
  void parse_gradient_flow_measure(Yp::inspect_node &in,
                                   const std::vector<std::string> &inner_tree,
                                   gp::measure_gradient_flow &mgfparams) {
    const std::vector<std::string> state0 = in.get_InnerTree();
    in.dig_deeper(inner_tree); // entering the glueball node
    YAML::Node nd = in.get_outer_node();

    mgfparams.measure_it = true;
    in.read_verb<double>(mgfparams.epsilon, {"epsilon"});
    in.read_verb<double>(mgfparams.tmax, {"tmax"});
    in.read_opt_verb<double>(mgfparams.tstart, {"tstart"});
    in.read_verb<bool>(mgfparams.save_conf, {"save_conf"});

    in.set_InnerTree(state0); // reset to previous state
  }

  /**
   * @brief parsing the online measurement block of the YAML input file
   *
   * @param in inspection node (full tree)
   * @param inner_tree path to the given branch of the tree as a std::vector of names of
   * branches
   * @param mparams reference to the measurements parameters
   * @param type "online" or "offline"
   */
  void parse_omeas(Yp::inspect_node &in,
                   const std::vector<std::string> &inner_tree,
                   gp::measure &mparams) {
    const std::vector<std::string> state0 = in.get_InnerTree();
    in.dig_deeper(inner_tree); // entering the glueball node
    YAML::Node nd = in.get_outer_node();

    mparams.do_mcmc = !bool(nd["offline"]); // true if doing offline measurements
    if (!mparams.do_mcmc) { // if offline, need these parameters
      in.read_verb<std::string>(mparams.conf_dir, {"offline", "conf_dir"});
      in.read_opt_verb<std::string>(mparams.conf_basename, {"offline", "conf_basename"});
      in.read_opt_verb<bool>(mparams.lenghty_conf_name, {"offline", "lenghty_conf_name"});
      in.read_opt_verb<size_t>(mparams.beta_str_width, {"offline", "beta_str_width"});
      validate_beta_str_width(mparams.beta_str_width);
    }

    mparams.res_dir = mparams.conf_dir; // default
    in.read_opt_verb<std::string>(mparams.res_dir, {"res_dir"});

    in.read_opt_verb<size_t>(mparams.verbosity, {"verbosity"});

    in.read_opt_verb<size_t>(mparams.icounter, {"icounter"});
    in.read_opt_verb<size_t>(mparams.nstep, {"nstep"});
    in.read_opt_verb<size_t>(mparams.n_meas, {"n_meas"});

    if (nd["pion_staggered"]) {
      mparams.pion_staggered = true;
      in.read_verb<double>(mparams.m0, {"pion_staggered", "mass"});
    }

    in.read_opt_verb<bool>(mparams.Wloop, {"Wloop"});

    // optional parameters for potentials
    if (nd["potential"]) {
      in.read_opt_verb<bool>(mparams.potentialplanar, {"potential", "potentialplanar"});
      in.read_opt_verb<bool>(mparams.potentialnonplanar,
                             {"potential", "potentialnonplanar"});
      in.read_opt_verb<bool>(mparams.append, {"potential", "append"});
      in.read_opt_verb<bool>(mparams.smear_spatial_only,
                             {"potential", "smear_spatial_only"});
      in.read_opt_verb<bool>(mparams.smear_temporal_only,
                             {"potential", "smear_temporal_only"});
      in.read_opt_verb<size_t>(mparams.n_apesmear, {"potential", "n_apesmear"});
      in.read_opt_verb<double>(mparams.alpha, {"potential", "alpha"});
      in.read_opt_verb<double>(mparams.sizeWloops, {"potential", "sizeWloops"});
    }
    if (nd["glueball"]) {
      parse_glueball_measure(in, {"glueball"}, mparams.glueball);
    }

    if (nd["gradient_flow"]) {
      parse_gradient_flow_measure(in, {"gradient_flow"}, mparams.gradient_flow);
    }

    in.set_InnerTree(state0); // reset to previous state
  }

  /**
   * @brief parsing the `metropolis` block of the YAML input file
   *
   * @param in inspection node (full tree)
   * @param inner_tree path to the given branch of the tree
   * @param mcparams reference to the hmc parameters
   */
  void parse_metropolis(Yp::inspect_node &in,
                        const std::vector<std::string> &inner_tree,
                        gp::metropolis &mcparams) {
    const std::vector<std::string> state0 = in.get_InnerTree();
    in.dig_deeper(inner_tree); // entering the glueball node
    YAML::Node nd = in.get_outer_node();

    in.read_opt_verb<bool>(mcparams.do_mcmc, {"do_mcmc"});
    in.read_opt_verb<size_t>(mcparams.n_meas, {"n_meas"});
    in.read_opt_verb<size_t>(mcparams.N_save, {"N_save"});
    in.read_opt_verb<size_t>(mcparams.seed, {"seed"});

    in.read_opt_verb<std::string>(mcparams.conf_dir, {"conf_dir"});
    in.read_opt_verb<std::string>(mcparams.conf_basename, {"conf_basename"});
    in.read_opt_verb<bool>(mcparams.lenghty_conf_name, {"lenghty_conf_name"});
    in.read_opt_verb<size_t>(mcparams.beta_str_width, {"beta_str_width"});
    validate_beta_str_width(mcparams.beta_str_width);

    heat_restart_incompatible(nd);

    in.read_opt_verb<bool>(mcparams.restart, {"restart"});
    in.read_opt_verb<bool>(mcparams.heat, {"heat"});

    in.read_verb<double>(mcparams.delta, {"delta"});
    in.read_opt_verb<size_t>(mcparams.N_hit, {"N_hit"});
    validate_N_hit(mcparams.N_hit);

    in.set_InnerTree(state0); // reset to previous state
    return;
  }

  /**
   * @brief parsing the hmc block of the YAML input file
   *
   * @param in inspection node (full tree)
   * @param inner_tree path to the given branch of the tree as a std::vector of names of
   * branches
   * @param hparams reference to the hmc parameters
   */
  void parse_hmc(Yp::inspect_node &in,
                 const std::vector<std::string> &inner_tree,
                 gp::hmc &hparams) {
    const std::vector<std::string> state0 = in.get_InnerTree();
    in.dig_deeper(inner_tree); // entering the glueball node
    YAML::Node nd = in.get_outer_node();

    in.read_verb<size_t>(hparams.N_save, {"n_save"});
    in.read_verb<size_t>(hparams.n_meas, {"n_meas"});

    in.read_opt_verb<bool>(hparams.do_mcmc, {"do_mcmc"});

    heat_restart_incompatible(nd);

    in.read_opt_verb<bool>(hparams.restart, {"restart"});
    in.read_opt_verb<bool>(hparams.heat, {"heat"});

    in.read_opt_verb<size_t>(hparams.seed, {"seed"});
    in.read_opt_verb<std::string>(hparams.configfilename, {"configname"});
    in.read_opt_verb<std::string>(hparams.conf_dir, {"conf_dir"});
    in.read_opt_verb<std::string>(hparams.conf_basename, {"conf_basename"});
    in.read_opt_verb<bool>(hparams.lenghty_conf_name, {"lenghty_conf_name"});

    in.read_opt_verb<size_t>(hparams.beta_str_width, {"beta_str_width"});
    validate_beta_str_width(hparams.beta_str_width);

    in.set_InnerTree(state0); // reset to previous state
    return;
  }

  void parse_integrator(Yp::inspect_node &in,
                        const std::vector<std::string> &inner_tree,
                        gp::hmc &hparams) {
    const std::vector<std::string> state0 = in.get_InnerTree();
    in.dig_deeper(inner_tree); // entering the glueball node
    YAML::Node nd = in.get_outer_node();

    in.read_opt_verb<size_t>(hparams.N_rev, {"N_rev"});
    in.read_opt_verb<size_t>(hparams.n_steps, {"n_steps"});
    in.read_opt_verb<double>(hparams.tau, {"tau"});
    in.read_opt_verb<size_t>(hparams.exponent, {"exponent"});
    in.read_opt_verb<std::string>(hparams.integrator, {"name"});

    in.set_InnerTree(state0); // reset to previous state
  }

  /**
   * @brief parse monomials and operators
   *
   * @param in
   * @param inner_tree
   * @param pparams reference to the physics parameters
   * @param sparams reference to the parameters specific to the program
   */
  template <class S>
  void parse_action(Yp::inspect_node &in,
                    const std::vector<std::string> &inner_tree,
                    gp::physics &pparams,
                    S &sparams) {
    const std::vector<std::string> state0 = in.get_InnerTree();
    in.dig_deeper(inner_tree); // entering the glueball node
    YAML::Node nd = in.get_outer_node();

    if (nd["monomials"]) {
      if (nd["monomials"]["gauge"]) {
        pparams.include_gauge = true;
        in.read_verb<double>(pparams.beta, {"monomials", "gauge", "beta"});
        if (nd["monomials"]["gauge"]["anisotropic"]) {
          pparams.anisotropic = true;
          in.read_opt_verb<double>(pparams.xi,
                                   {"monomials", "gauge", "anisotropic", "xi"});
        }
      }
      if (nd["operators"]) {
        // note: initializing the D*Ddag monomial makes sense only with a valid Dirac
        // operator
        if (nd["operators"]["staggered"] && nd["monomials"]["staggered_det_DDdag"]) {
          pparams.include_staggered_fermions = true;

          in.read_verb<double>(pparams.m0, {"operators", "staggered", "mass"});

          in.read_opt_verb<std::string>(sparams.solver,
                                        {"monomials", "staggered_det_DDdag", "solver"});
          in.read_opt_verb<double>(sparams.tolerance_cg,
                                   {"monomials", "staggered_det_DDdag", "tolerance_cg"});
          in.read_opt_verb<size_t>(
            sparams.solver_verbosity,
            {"monomials", "staggered_det_DDdag", "solver_verbosity"});
          in.read_opt_verb<size_t>(sparams.seed_pf,
                                   {"monomials", "staggered_det_DDdag", "seed_pf"});
        }
      }
    }

    in.set_InnerTree(state0); // reset to previous state
  }

  namespace hmc {
    void parse_input_file(const YAML::Node &nd, gp::physics &pparams, gp::hmc &hparams) {
      Yp::inspect_node in(nd);

      parse_geometry(in, pparams);
      parse_action<gp::hmc>(in, {}, pparams, hparams);

      parse_hmc(in, {"hmc"}, hparams); // hmc-u1 parameters
      parse_integrator(in, {"integrator"}, hparams); // integrator parameters

      // online measurements
      if (nd["omeas"]) {
        hparams.do_omeas = true;
        hparams.omeas.conf_dir = hparams.conf_dir;
        parse_omeas(in, {"omeas"}, hparams.omeas);
      }

      in.finalize();
      return;
    }

  } // namespace hmc

  namespace measure {

    void
    parse_input_file(const YAML::Node &nd, gp::physics &pparams, gp::measure &mparams) {
      Yp::inspect_node in(nd);

      parse_geometry(in, pparams);

      // beta, xi value from the gauge action
      in.read_verb<double>(pparams.beta, {"monomials", "gauge", "beta"});
      if (nd["monomials"]["gauge"]["anisotropic"]) {
        pparams.anisotropic = true;
        in.read_opt_verb<double>(pparams.xi, {"monomials", "gauge", "anisotropic", "xi"});
      }

      if (nd["omeas"]) {
        parse_omeas(in, {"omeas"}, mparams);
      }

      in.finalize();
      return;
    }

  } // namespace measure

  namespace metropolis {

    void parse_input_file(const YAML::Node &nd,
                          gp::physics &pparams,
                          gp::metropolis &mcparams) {
      Yp::inspect_node in(nd);

      parse_geometry(in, pparams);
      parse_action<gp::metropolis>(in, {}, pparams, mcparams);
      parse_metropolis(in, {"metropolis"}, mcparams);

      if (nd["omeas"]) {
        mcparams.do_omeas = true;
        mcparams.omeas.conf_dir = mcparams.conf_dir;
        parse_omeas(in, {"omeas"}, mcparams.omeas);
      }

      in.finalize();

      return;
    }

  } // namespace metropolis

} // namespace input_file_parsing
