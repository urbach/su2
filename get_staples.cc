#include"su2.hh"
#include"gaugeconfig.hh"
#include"get_staples.hh"
#include<vector>

//
// U^staples = \sum_{nu != mu} U_nu(x+mu) U_mu(x+nu)^dagger + U_nu(x)^dagger
//             + sum__{nu != mu} U_nu(x+mu-nu)^dagger U_mu(x-nu)^dagger U_nu(x-nu)
// 

su2 get_staples(gaugeconfig &U, vector<size_t> const x, const size_t mu) {
  su2 K(0., 0.);
  vector<size_t> x1 = x, x2 = x;
  x1[mu] += 1;
  for(size_t nu = 0; nu < U.getndims(); nu++) {
    if(nu != mu) {
      x2[nu]++;
      K += U(x1, nu) * U(x2, mu).dagger() * U(x, nu).dagger();
      x2[nu]--;
    }
  }
  for(size_t nu = 0; nu < U.getndims(); nu++) {
    if(nu != mu) {
      x1[nu]--;
      x2[nu]--;
      K += U(x1, nu).dagger() * U(x2, mu).dagger() * U(x2, nu);
      x2[nu]++;
      x1[nu]++;
    }
  }
  return(K);
}
