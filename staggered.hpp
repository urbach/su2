// staggered.hpp
/**
 * Staggered fermions 
 */

#pragma once

#include <array>
#include <vector>
#include <complex>
#include <random>

#include "adjointfield.hh"
#include "gauge_energy.hh"
#include "gaugeconfig.hh"
#include "get_staples.hh"
#include "hamiltonian_field.hh"
#include "monomial.hh"
#include "su2.hh"
#include "u1.hh"
#include "include/geometry.hh"

#include "cg.hpp"


namespace staggered {

// eta_{\mu}(x) as in eq. (16) of https://arxiv.org/pdf/2112.14640.pdf
// or eq. (8) of https://www.sciencedirect.com/science/article/pii/0550321389903246
double eta(const std::vector<size_t>& x, const size_t& mu){
  double s = 0;
  for(int nu=0; nu<mu; nu++) { s += x[nu]; }
  return std::pow(-1.0, s);
}

// index of the lattice point given the dimensions
size_t txyz_to_index(
  const size_t& t, const size_t& x, const size_t& y, const size_t& z,
  const size_t& Lt, const size_t& Lx, const size_t& Ly, const size_t& Lz
  ) 
{
  const geometry g(Lx, Ly, Lz, Lt); // note the order
  return g.getIndex(t, x, y, z);
}

// overload : x={x0,x1,x2,x3}, dims={Lt,Lx,Ly,Lz}
size_t txyz_to_index(const std::vector<size_t>& x, const std::vector<size_t>& dims) 
{
  return txyz_to_index(x[0], x[1], x[2], x[3], dims[0], dims[1], dims[2], dims[3]);
}

// vector of staggered "spinors" (no Dirac structure) for all the points of the lattice
template<class Float, class Type>
class spinor_lat_4d {
private:
std::vector<Type> Psi;

public:

spinor_lat_4d(){}
~spinor_lat_4d(){}

Type& operator () (const size_t& i){ return Psi[i]; }
Type operator () (const size_t& i) const { return Psi[i]; }


Type& operator () (const std::vector<size_t>& x, const std::vector<size_t>& dims){
  const size_t i = txyz_to_index(x, dims);
  return (*this)[i];
}

Type operator () (const std::vector<size_t>& x, const std::vector<size_t>& dims) const {
  const size_t i = txyz_to_index(x, dims);
  return (*this)[i];
}

Type& operator [] (const size_t& i){ return Psi[i];}
Type operator [] (const size_t& i) const { return Psi[i]; }

spinor_lat_4d(const size_t& n){ Psi.resize(n); }

spinor_lat_4d(const size_t& n, const Type& val){
  std::vector<Type> v(n, val);
  (*this).Psi = v;  
}

size_t size() const { return Psi.size();}
void resize(const size_t& n, const Type& val = (Type) 0.0){ Psi.resize(n, val); }

spinor_lat_4d<Float, Type> operator /(const Type& lambda){
  const int N = (*this).size();
  spinor_lat_4d<Float, Type> phi(N);
  for (size_t i = 0; i < N; i++) { phi[i] = Psi[i]/lambda; }
  return phi;
}


void operator +=(const spinor_lat_4d<Float, Type>& psi){ (*this) = (*this) + psi; }

Float norm_squared() const{ 
  return complex_dot_product((*this),(*this)).real(); // A^{\dagger} * A is real
}

Float norm() const { return sqrt(this->norm_squared()); }

}; // class spinor_lat_4d


// a + lambda*b
template<class Float, class Type, class Type_lambda>
spinor_lat_4d<Float, Type> a_plus_lambda_b(const spinor_lat_4d<Float, Type>& a, const Type_lambda& lambda, const spinor_lat_4d<Float, Type>& b) 
{
  const int N = a.size();
  spinor_lat_4d<Float, Type> c(N);
  for (size_t i = 0; i < N; i++) { c[i] = a[i] + lambda*b[i]; }
  return c;
}


template<class Float, class Type>
spinor_lat_4d<Float, Type> operator +(const spinor_lat_4d<Float, Type>& a, const spinor_lat_4d<Float, Type>& b){
  return a_plus_lambda_b(a, 1.0, b);
}

template<class Float, class Type>
spinor_lat_4d<Float, Type> operator -(const spinor_lat_4d<Float, Type>& a, const spinor_lat_4d<Float, Type>& b) {
  return a_plus_lambda_b(a, -1.0, b);
}

// change of sign : psi --> -psi
template<class Float, class Type>
spinor_lat_4d<Float, Type> operator -(const spinor_lat_4d<Float, Type>& psi){
  const spinor_lat_4d<Float, Type> v( psi.size() );
  return (v-psi);
}


template<class Float, class Type >
spinor_lat_4d<Float, Type> gaussian_spinor_normalized(const size_t& n, const Float& avr, const Float& sigma, const size_t& seed)
{
  std::normal_distribution<Float> dis{avr, sigma};
  std::mt19937 gen_re(seed), gen_im(seed+1);
  spinor_lat_4d<Float, Type> psi_gauss(n);
  Type x, y;
  Type norm2 = 0.0;
  for (size_t i = 0; i < n; i++) {   // lattice points
    x = (dis(gen_re), dis(gen_im));
    psi_gauss[i] = x;
    norm2 += (conj(x)*x).real(); // x^{\dagger}*x is real
  }
  return psi_gauss / sqrt(norm2);
}



// \sum_{i} A_i^{\dagger}*B_i
template<class Float, class Type>
Type complex_dot_product(const spinor_lat_4d<Float, Type>&A, const spinor_lat_4d<Float, Type>&B)
{
  const int N = A.size();
  Type sum = 0.0;
  
  // #pragma omp parallel for
  for (size_t i = 0; i < N; i++) { 
    sum += conj(A[i]) * B[i]; 
  }

  return sum;
}

template<class Float, class Type>
Type operator *(const spinor_lat_4d<Float, Type>&A, const spinor_lat_4d<Float, Type>&B)
{ return complex_dot_product(A, B); }


template<class Float, class Type>
spinor_lat_4d<Float, Type> operator *(const Type& lambda, const spinor_lat_4d<Float, Type>& psi)
{
  const spinor_lat_4d<Float, Type> v(psi.size());
  return a_plus_lambda_b(v, lambda, psi);
}

template<class Float, class Type>
Float norm(const spinor_lat_4d<Float, Type>& psi){ return psi.norm(); }


// class for D^{\dagger}*D
template<class Float, class Type, class Group>
class DdagD_matrix_lat{
  public:
  gaugeconfig<Group>* U;  
  Float m;

  DdagD_matrix_lat(){}
  ~DdagD_matrix_lat(){}

  DdagD_matrix_lat(gaugeconfig<Group>* _U, const Float& _m)
  { 
    U = _U;
    m = _m;
  }

  // rows = cols = N = number of lattice points
  size_t rows() const { return U->getVolume(); }
  size_t cols() const { return this->rows(); }


  spinor_lat_4d<Float, Type> inv(const spinor_lat_4d<Float, Type> & psi, const Float& tol, const size_t& verb, const size_t& seed) const {
  typedef spinor_lat_4d<Float, Type> LAvector;
  typedef DdagD_matrix_lat<Float, Type, Group> LAmatrix;
  cg::LinearCG<Float, Type, LAmatrix, LAvector>  LCG((*this), psi);
  const size_t N = psi.size();
  const LAvector phi0 = staggered::gaussian_spinor_normalized<Float, Complex>(N, 0.0, 1.0, seed);
 
  if(verb>1){ std::cout << "Calling the CG solver.\n"; }

  LCG.solve(phi0, tol, verb);
  return LCG.get_solution();
  }

};

template<class Float, class Type, class Group>
spinor_lat_4d<Float, Type> operator *(
  const DdagD_matrix_lat<Float, Type, Group>& M, 
  const spinor_lat_4d<Float, Type>& psi)
{
  gaugeconfig<Group>* U = M.U;
  const Float m = M.m;

  const size_t Lt = U->getLt(), Lx = U->getLx(), Ly = U->getLy(), Lz = U->getLz();
  const std::vector<size_t> dims = {Lt, Lx, Ly, Lz}; // vector of dimensions
  const size_t nd = U->getndims();

  const int N = psi.size();
  spinor_lat_4d<Float, Type> phi(N);
#pragma omp parallel for
  for(size_t x0 = 0; x0 < Lt; x0++) {
    for(size_t x1 = 0; x1 < Lx; x1++) {
      for(size_t x2 = 0; x2 < Ly; x2++) {
        for(size_t x3 = 0; x3 < Lz; x3++) {
          const std::vector<size_t> x = {x0, x1, x2, x3};
          std::vector<size_t> xm = x; // will be x+mu
          std::vector<size_t> xp = x; // sill be x-mu
          std::vector<size_t> xpm = x; // will be x + mu - nu
          std::vector<size_t> xmm = x; // will be x - mu - nu
          std::vector<size_t> xmp = x; // will be  x - mu + nu
          std::vector<size_t> xpp = x; // will be  x + mu + nu
          for(size_t mu = 0; mu < nd; mu++) {
            xp[mu] += 1; // x + mu
            xm[mu] -= 1; // x - mu
            xpp[mu] += 1; // see later in the loop
            xpm[mu] += 1; // see later in the loop
            xmp[mu] -= 1; // see later in the loop
            xmm[mu] -= 1; // see later in the loop
            //
            const Float fact_mu   = (1.0/2.0) * eta(x, mu);
            for(size_t nu = 0; nu < nd; nu++) {
              xpp[nu] += 1; // x+mu+nu
              xpm[nu] -= 1; // x+mu-nu
              xmp[nu] += 1; // x-mu+nu
              xmm[nu] -= 1; // x-mu-nu

              const Float fact_nu_p = (1.0/2.0) * eta(xp, nu);
              const Float fact_nu_m = (1.0/2.0) * eta(xm, nu);

              const Float fact1 = fact_mu * fact_nu_p;
              const Float fact2 = fact_mu * fact_nu_m;

              phi(x, dims) += 
                fact1 * (*U)(x, mu).dagger() * (*U)(xp, mu) * psi(xpp, dims);
              phi(x, dims) -= 
                fact1 * (*U)(x, mu).dagger() * (*U)(xpm, mu).dagger() * psi(xmm, dims);
              phi(x, dims) -= 
                fact2 * (*U)(xm, mu) * (*U)(xm, mu) * psi(xmm, dims);
              phi(x, dims) += fact2 * (*U)(xm, mu) * (*U)(xmm, mu).dagger() * psi(xmm, dims);

              xpp[nu] -= 1; // x+mu again
              xpm[nu] += 1; // x+mu again
              xmp[nu] -= 1; // x-mu again
              xmm[nu] += 1; // x-mu again
            }
            //
            phi(x, dims) += m * fact_mu * ( (*U)(x, mu) + (*U)(x, mu).dagger() ) * psi(xp, dims);
            phi(x, dims) -= m * fact_mu * ( (*U)(xm, mu).dagger() + (*U)(xm, mu) ) * psi(xm, dims);

            xm[mu] += 1; // =x again
            xp[mu] -= 1; // =x again
            xpp[mu] -= 1; // =x again
            xpm[mu] -= 1; // =x again
            xmp[mu] += 1; // =x again
            xmm[mu] += 1; // =x again
          }
          //
          phi(x, dims) = std::pow(m, 2.0) * psi(x, dims);
        }
      }
    }
  }
  //
  return phi;  
}




// returns the reslut of D*psi, where D is the Dirac operator
template<class Float, class Type, class Group>
spinor_lat_4d<Float, Type> apply_D(gaugeconfig<Group>* U, const Float& m, const spinor_lat_4d<Float, Type>& psi)
{
  const size_t Lt = U->getLt(), Lx = U->getLx(), Ly = U->getLy(), Lz = U->getLz();
  const std::vector<size_t> dims = {Lt, Lx, Ly, Lz}; // vector of dimensions
  const size_t nd = U->getndims();

  const int N = psi.size();
  spinor_lat_4d<Float, Type> phi(N);
#pragma omp parallel for
  for(size_t x0 = 0; x0 < Lt; x0++) {
    for(size_t x1 = 0; x1 < Lx; x1++) {
      for(size_t x2 = 0; x2 < Ly; x2++) {
        for(size_t x3 = 0; x3 < Lz; x3++) {
          const std::vector<size_t> x = {x0, x1, x2, x3};
          std::vector<size_t> xm = x, xp = x;
          for(size_t mu = 0; mu < nd; mu++) {
            const Float eta_x_mu = eta(x, mu);
            xm[mu] -= 1; // x - mu
            xp[mu] += 1; // x + mu
            
            phi(x, dims) += (1.0/2.0) * eta_x_mu * (*U)(x, mu) * psi(xp, dims);
            phi(x, dims) -= (1.0/2.0) * eta_x_mu * (*U)(x, mu).dagger() * psi(xm, dims);

            xm[mu] += 1; // =x again
            xp[mu] -= 1; // =x again
          }
          phi(x, dims) += m * psi(x, dims);
        }
      }
    }
  }
  return phi;  
} 


// returns the reslut of D^{\dagger}*psi, where D is the Dirac operator
template<class Float, class Type, class Group>
spinor_lat_4d<Float, Type> apply_Ddag(gaugeconfig<Group>* U, const Float& m, const spinor_lat_4d<Float, Type>& psi)
{
  const size_t Lt = U->getLt(), Lx = U->getLx(), Ly = U->getLy(), Lz = U->getLz();
  const std::vector<size_t> dims = {Lt, Lx, Ly, Lz}; // vector of dimensions
  const size_t nd = U->getndims();

  const int N = psi.size();
  spinor_lat_4d<Float, Type> phi(N);
#pragma omp parallel for
  for(size_t x0 = 0; x0 < Lt; x0++) {
    for(size_t x1 = 0; x1 < Lx; x1++) {
      for(size_t x2 = 0; x2 < Ly; x2++) {
        for(size_t x3 = 0; x3 < Lz; x3++) {
          const std::vector<size_t> x = {x0, x1, x2, x3};
          std::vector<size_t> xm = x, xp = x;
          for(size_t mu = 0; mu < nd; mu++) {
            const Float eta_x_mu = eta(x, mu);
            xm[mu] -= 1; // x - mu
            xp[mu] += 1; // x + mu
            
            phi(x, dims) += (1.0/2.0) * eta_x_mu * (*U)(x, mu).dagger() * psi(xp, dims);
            phi(x, dims) -= (1.0/2.0) * eta_x_mu * (*U)(x, mu) * psi(xm, dims);

            xm[mu] += 1; // =x again
            xp[mu] -= 1; // =x again
          }
          phi(x, dims) += m * psi(x, dims);
        }
      }
    }
  }
  return phi;  
} 


// returns the result of  dD*psi, where dD is the derivative of the Dirac operator for a U(1) theory 
// as in eq. 8.36 of Gattringer&Lang
template<class Float, class Type, class Group>
spinor_lat_4d<Float, Type> apply_der_D(const std::vector<size_t>& x, const size_t& mu, gaugeconfig<Group>* U, const Float& m, const spinor_lat_4d<Float, Type>& psi)
{
  const size_t Lt = U->getLt(), Lx = U->getLx(), Ly = U->getLy(), Lz = U->getLz();
  const std::vector<size_t> dims = {Lt, Lx, Ly, Lz}; // vector of dimensions

  spinor_lat_4d<Float, Type> phi(psi.size());

  std::vector<size_t> xm = x, xp = x;
  xm[mu]--; // x - mu
  xp[mu]++; // x + mu

  const Float eta_x_mu = eta(x, mu);
  const Float eta_xm_mu = eta(xm, mu);
  const std::complex<Float> i(0.0, 1.0);
  phi(x, dims) = (1.0/2.0) * eta_x_mu * (+i) * (*U)(x, mu) * psi(xp, dims);
  phi(xm, dims) = -(1.0/2.0) * eta_xm_mu * (-i) * (*U)(x, mu).dagger() * psi(xm, dims);
  
  return phi;
}


// returns the reslut of  dD^{\dagger}*psi, where dD^{\dagger} is the derivative of the Dirac operator (daggered) for a U(1) theory 
// analogously to in eq. 8.36 of Gattringer&Lang
template<class Float, class Type, class Group>
spinor_lat_4d<Float, Type> apply_der_Ddag(const std::vector<size_t>& x, const size_t& mu, gaugeconfig<Group>* U, const Float& m, const spinor_lat_4d<Float, Type>& psi)
{
  const size_t Lt = U->getLt(), Lx = U->getLx(), Ly = U->getLy(), Lz = U->getLz();
  const std::vector<size_t> dims = {Lt, Lx, Ly, Lz}; // vector of dimensions

  spinor_lat_4d<Float, Type> phi(psi.size());

  std::vector<size_t> xm = x, xp = x;
  xm[mu]--; // x - mu
  xp[mu]++; // x + mu

  const Float eta_x_mu = eta(x, mu);
  const Float eta_xm_mu = eta(xm, mu);
  const std::complex<Float> i(0.0, 1.0);
  phi(x, dims)  = (1.0/2.0) * eta_x_mu * (-i) * (*U)(x, mu).dagger() * psi(xp, dims);
  phi(xm, dims) = -(1.0/2.0) * eta_xm_mu * (+i) * (*U)(x, mu) * psi(xm, dims);
  
  return phi;
} 



} // namespace staggered
