/**
 * @file heatbath_overrelaxation.hpp
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class for heatbath_overrelaxation Algorithm
 * @version 0.1
 * @date 2022-09-01
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "base_program.hpp"
#include "heatbath.hh"
#include "overrelaxation.hh"
#include "parse_input_file.hh"

template <class Group>
class heatbath_overrelaxation_algo
  : public base_program<Group, gp::heatbath_overrelaxation> {
private:
  std::vector<double> rate = {0.0, 0.0, 0.0, 0.0};

public:
  heatbath_overrelaxation_algo() { (*this).algo_name = "heatbath_overrelaxation"; }
  ~heatbath_overrelaxation_algo() {}

  void print_program_info() const {
    std::cout << "## Heatbath+Overrelaxation Algorithm\n";
  }

  void parse_input_file(const YAML::Node &nd) {
    namespace in_heatbath_overrelaxation = input_file_parsing::heatbath_overrelaxation;
    in_heatbath_overrelaxation::parse_input_file(nd, (*this).pparams, (*this).sparams);
    (*this).omeas = (*this).sparams.omeas;
    (*this).conf_path_basename =
      io::get_conf_path_basename((*this).pparams, (*this).sparams);
  }

  void set_omp_threads() {
#ifdef _USE_OMP_
    /**
     * the parallelisation of the sweep-function first iterates over all odd points in t
     * and then over all even points because the nearest neighbours must not change
     * during the updates, this is not possible for an uneven number of points in T
     * */
    if ((*this).pparams.Lt % 2 != 0) {
      std::cerr << "For parallel computing an even number of points in T is needed!"
                << std::endl;
      omp_set_num_threads(1);
      std::cerr << "Continuing with one thread." << std::endl;
    }
    // set things up for parallel computing in sweep
    (*this).threads = omp_get_max_threads();
#else
    (*this).threads = 1;
#endif
    std::cout << "threads " << (*this).threads << std::endl;
  }

  /**
   * @brief do the i-th sweep of the heatbath_overrelaxation algorithm
   *
   * @param i trajectory index
   */
  void do_heatbath(const size_t &i, const size_t &i_hb) {
    if ((*this).sparams.do_mcmc) {
      const int n_threads = (*this).threads;
      std::vector<std::mt19937> engines(n_threads);
      for (size_t i_engine = 0; i_engine < n_threads; i_engine++) {
        engines[i_engine].seed(((*this).sparams.seed + i * n_threads + i_engine) *
                                 (*this).sparams.n_heatbath +
                               i_hb);
      }

      (*this).rate += heatbath((*this).U, engines, (*this).pparams.beta,
                               (*this).pparams.xi, (*this).pparams.anisotropic);
    }
  }

  void do_overrelaxation() {
    overrelaxation((*this).U, (*this).pparams.beta, (*this).pparams.xi,
                   (*this).pparams.anisotropic);
  }

  // save acceptance rates to additional file to keep track of measurements
  void save_acceptance_rates() {
    if ((*this).sparams.do_mcmc) {
      std::cout << "## Acceptanced links "
                << rate[0] / double((*this).sparams.n_meas * (*this).sparams.n_heatbath)
                << " accepted temporal links "
                << rate[1] / double((*this).sparams.n_meas * (*this).sparams.n_heatbath)
                << " acceptance rate "
                << rate[2] / double((*this).sparams.n_meas * (*this).sparams.n_heatbath)
                << " temporal acceptance rate "
                << rate[3] / double((*this).sparams.n_meas * (*this).sparams.n_heatbath)
                << std::endl;
      (*this).acceptancerates.open((*this).sparams.conf_dir +
                                     "/acceptancerates-heatbath_overrelaxation.data",
                                   std::ios::app);
      (*this).acceptancerates
        << rate[0] / double((*this).sparams.n_meas * (*this).sparams.n_heatbath) << " "
        << rate[1] / double((*this).sparams.n_meas * (*this).sparams.n_heatbath) << " "
        << rate[2] / double((*this).sparams.n_meas * (*this).sparams.n_heatbath) << " "
        << rate[3] / double((*this).sparams.n_meas * (*this).sparams.n_heatbath) << " "
        << (*this).pparams.beta << " " << (*this).pparams.Lx << " " << (*this).pparams.Lt
        << " " << (*this).pparams.xi << " " << (*this).sparams.heat << " "
        << (*this).threads << " " << (*this).sparams.n_meas << " " << (*this).sparams.seed
        << " " << std::endl;
      (*this).acceptancerates.close();
    }
    return;
  }

  void run(const YAML::Node &nd) {
    this->pre_run(nd);
    this->init_gauge_conf_mcmc();
    this->open_output_data();
    this->set_omp_threads();
    if ((*this).sparams.do_omeas) {
      this->set_potential_filenames();
    }

    if ((*this).g_icounter == 0) {
      // header: column names in the output
      std::string head_str = io::get_header_1(" ");
      std::cout << head_str;
      (*this).os << head_str;
    }

    size_t i_min = (*this).g_icounter;
    size_t i_max = (*this).sparams.n_meas + (*this).g_icounter;

    /**
     * do measurements:
     * heatbath_overrelaxation updates of every link in the lattice
     * calculate plaquette, spacial plaquette, energy density with and without cloverdef
     * and write to stdout and output-file save every nave configuration
     * */
    const int seed = (*this).sparams.seed;
    const int n_threads = (*this).threads;
    std::vector<std::mt19937> engines(n_threads);
    for (size_t i = i_min; i < i_max; i++) {
      if ((*this).sparams.do_mcmc) {
        this->output_line(i);

        const int d1 = i;
        const int d2 = n_threads;
        const int d3 = (*this).sparams.n_heatbath;
        for (size_t i_hb = 0; i_hb < (*this).sparams.n_heatbath; i_hb++) {
          for (size_t i_engine = 0; i_engine < n_threads; i_engine++) {
            engines[i_engine].seed(seed + i * d2 * d3 + i_engine * d3 + i_hb);
          }
          this->do_heatbath(i, i_hb);
        }

        for (size_t over = 0; over < (*this).sparams.n_overrelax; over++) {
          this->do_overrelaxation();
        }

        if (i > 0 && (i % (*this).sparams.N_save) == 0) {
          std::ostringstream oss_i;
          oss_i << (*this).conf_path_basename << "." << i << std::ends;
          ((*this).U).save(oss_i.str());
        }
      }

      bool b1 = (*this).sparams.do_omeas;
      bool b2 = i != 0;
      bool b3 = (i % (*this).sparams.N_save) == 0;
      bool do_omeas = (b1 && b2 && b3);
      this->after_MCMC_step(i, do_omeas);
    }
    this->save_acceptance_rates();
    this->save_final_conf();
  }
};
