/**
 * @file gaugeconfig.hh
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief class for U(1) or SU(2) gauge configuration
 * @version 0.1
 * @date 2022-10-03
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include "geometry.hh"
#include "random_element.hh"
#include "su2.hh"
#include "u1.hh"

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

// test function to see if the template has been specified correctly
template <class T> std::string gauge_group_name();
template <> std::string gauge_group_name<_u1>() {
  return "u1";
}
template <> std::string gauge_group_name<_su2>() {
  return "su2";
}

using std::vector;

template <class T> class gaugeconfig {
  template <class iT> using nd_max_arr = spacetime_lattice::nd_max_arr<iT>;

public:
  using value_type = T;

  gaugeconfig(){};
  ~gaugeconfig(){};

  gaugeconfig(const size_t Lx,
              const size_t Ly,
              const size_t Lz,
              const size_t Lt,
              const size_t ndims = spacetime_lattice::nd_max,
              const double beta = 0)
    : Lx(Lx),
      Ly(Ly),
      Lz(Lz),
      Lt(Lt),
      volume(Lx * Ly * Lz * Lt),
      beta(beta),
      ndims(ndims) {
    data.resize(volume * ndims);
  }
  gaugeconfig(const gaugeconfig &U)
    : Lx(U.getLx()),
      Ly(U.getLy()),
      Lz(U.getLz()),
      Lt(U.getLt()),
      volume(U.getVolume()),
      beta(U.getBeta()),
      ndims(U.getndims()) {
    data.resize(volume * ndims);
#pragma omp parallel for
    for (size_t i = 0; i < getSize(); i++) {
      data[i] = U[i];
    }
  }

  size_t storage_size() const {
    return data.size() * sizeof(value_type);
  };
  std::vector<value_type> get_data() const {
    return data;
  }
  size_t getLx() const {
    return (Lx);
  }
  size_t getLy() const {
    return (Ly);
  }
  size_t getLz() const {
    return (Lz);
  }
  size_t getLt() const {
    return (Lt);
  }
  size_t getndims() const {
    return (ndims);
  }
  size_t getVolume() const {
    return (volume);
  }
  size_t getSize() const {
    return (volume * ndims);
  }
  double getBeta() const {
    return beta;
  }
  void setBeta(const double _beta) {
    beta = _beta;
  }
  int getNc() const {
    return (data[0].N_c);
  }
  void restoreSU() {
#pragma omp parallel for
    for (size_t i = 0; i < getSize(); i++) {
      data[i].restoreSU();
    }
  }

  void operator=(const gaugeconfig &U) {
    volume = U.getVolume();
    Lx = U.getLx();
    Ly = U.getLy();
    Lz = U.getLz();
    Lt = U.getLt();
    ndims = U.getndims();
    beta = U.getBeta();
    data.resize(U.getSize());
#pragma omp parallel for
    for (size_t i = 0; i < U.getSize(); i++) {
      data[i] = U[i];
    }
  }

  value_type &operator()(
    size_t const t, size_t const x, size_t const y, size_t const z, size_t const mu) {
    return data[getIndex(t, x, y, z, mu)];
  }

  value_type operator()(size_t const t,
                        size_t const x,
                        size_t const y,
                        size_t const z,
                        size_t const mu) const {
    return data[getIndex(t, x, y, z, mu)];
  }

  value_type &operator()(std::vector<size_t> const &coords, size_t const mu) {
    return data[getIndex(coords[0], coords[1], coords[2], coords[3], mu)];
  }

  value_type operator()(std::vector<size_t> const &coords, size_t const mu) const {
    return data[getIndex(coords[0], coords[1], coords[2], coords[3], mu)];
  }

  /**
   * access elements according to the convention of
   * https://link.springer.com/book/10.1007/978-3-642-01850-3, eq. (2.34)
   */
  template <class iT = int>
  value_type
  operator()(const nd_max_arr<iT> &x, const size_t &mu, const bool &bs = true) const {
    const value_type Ux_mu = data[getIndex(x[0], x[1], x[2], x[3], mu)];
    std::array<iT, spacetime_lattice::nd_max> xm = x;

    xm[mu]--;
    const value_type Udagx_mu = data[getIndex(xm[0], xm[1], xm[2], xm[3], mu)].dagger();

    if (bs) {
      return Ux_mu;
    } else {
      return Udagx_mu;
    }
    //    return ((double) bs) * Ux_mu + (1 - ((double)bs)) * Udagx_mu;
  }

  value_type &operator[](size_t const index) {
    return data[index];
  }

  value_type operator[](size_t const index) const {
    return data[index];
  }

  bool operator==(const gaugeconfig &U) const {
    typedef typename accum_type<T>::type accum;

    const size_t N = this->getSize();
    bool ret = true;
#pragma omp parallel for
    for (size_t i = 0; i < N; i++) {
      const bool ci = (accum(data[i]) == accum(U[i]));
      if (!ci) {
        ret = false;
      }
    }
    return ret;
  }

  void save(std::string const &path) const;
  int load(std::string const &path);

private:
  size_t Lx, Ly, Lz, Lt, volume, ndims;
  double beta;

  std::vector<value_type> data;

  size_t getIndex(const size_t t,
                  const size_t x,
                  const size_t y,
                  const size_t z,
                  const size_t _mu) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Lx) % Lx;
    size_t y2 = (y + Ly) % Ly;
    size_t y3 = (z + Lz) % Lz;
    size_t mu = (_mu + ndims) % ndims;
    return ((((y0 * Lx + y1) * Ly + y2) * Lz + y3) * ndims + mu);
  }
};

template <class T> void gaugeconfig<T>::save(std::string const &path) const {
  std::ofstream ofs(path, std::ios::out | std::ios::binary);
  ofs.write(reinterpret_cast<char const *>(data.data()), storage_size());
  ofs.close();
  return;
}

/**
 * @brief loading gauge configuration from binary file
 * ACHTUNG: The new gauge configuration must be initialized with the correct size!
 * The geometrical size of the lattices of the configuration stored in the file and the
 * one in which the former is loaded to must be the same!
 * @tparam T gauge group of the configuration
 * @param path location of the previously generated configuration
 * @return int 0=success and 1=failure to load
 */
template <class T> int gaugeconfig<T>::load(std::string const &path) {
  std::cout << "## Reading config from file " << path << std::endl;
  std::ifstream ifs(path, std::ios::in | std::ios::binary);
  if (ifs) {
    ifs.read(reinterpret_cast<char *>(data.data()), storage_size());
    ifs.close();
    return 0;
  } else {
    std::cerr << "Error: could not read file from " << path << std::endl;
  }
  return 1;
}

template <class T> void coldstart(gaugeconfig<T> &config) {
#pragma omp parallel for
  for (size_t i = 0; i < config.getSize(); i++) {
    config[i] = T(1., 0.);
  }
}

template <class T> void coldstart(gaugeconfig<_u1> &config) {
#pragma omp parallel for
  for (size_t i = 0; i < config.getSize(); i++) {
    config[i] = _u1(0.);
  }
}

/**
 * @brief Initialize the gauge configuration to either hot or cold start.
 * Each value of the configuration array is set equal to a random number distributed as
 * specified by `random_element()`.
 * @tparam T
 * @param config gauge configuration to be initialized
 * @param seed seed of the random number generator
 * @param _delta parameter of the random distribution
 */
template <class T>
void hotstart(gaugeconfig<T> &config, const int seed, const double _delta) {
  double delta = _delta;
  if (delta < 0.)
    delta = 0;
  if (delta > 1.)
    delta = 1.;
  std::mt19937 engine(seed);

  for (size_t i = 0; i < config.getSize(); i++) {
    random_element(config[i], engine, delta);
  }
}

/**
 * @brief copy for the homonymous function which takes `double _delta` as parameter.
 * @param hot : _delta = (double) hot
 */
template <class T>
void hotstart(gaugeconfig<T> &config, const int seed, const bool &hot) {
  const double delta = (double)hot;
  std::mt19937 engine(seed);
  for (size_t i = 0; i < config.getSize(); i++) {
    random_element(config[i], engine, delta);
  }
}
