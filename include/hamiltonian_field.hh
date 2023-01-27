#pragma once

#include "adjointfield.hh"
#include "gaugeconfig.hh"
#include <vector>

template <typename Float, class Group> struct hamiltonian_field {
  adjointfield<Float, Group> *momenta;
  gaugeconfig<Group> *U;

  // Our normalization convention of the generators is different than
  // https://arxiv.org/pdf/1006.4518.pdf this factor takes that into account when using
  // this class for the gradiend flow evolution of the gauge fields.
  double Fact_Nc_force_Luscher = 1.0 / 4.0; // default: non abelian case

  hamiltonian_field(adjointfield<Float, Group> &momenta, gaugeconfig<Group> &U)
    : momenta(&momenta), U(&U) {
    if (U.getNc() == 1) {
      Fact_Nc_force_Luscher = 1.0 / 2.0;
    }
  }
};
