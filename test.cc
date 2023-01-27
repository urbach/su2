#include"su2.hh"
#include"u1.hh"
#include"random_element.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hpp"
#include"random_gauge_trafo.hh"
#include"parse_commandline.hh"
#include"flat-energy_density.hh"

#include<iostream>
#include<vector>

int main() {
  size_t Lx = 8, Ly = 8, Lz = 8, Lt = 16;


  std::cout << "Tests of SU(2)" << std::endl << std::endl;;
  
  std::vector<su2> config;
  config.resize(Lx*Ly*Lz*Lt);
  
  // set all to 1
  for(auto i = config.begin(), end =config.end(); i < end; i++) {
    *i = su2(1., 0.);
  }

  su2 U(1., 0.);
  std::cout << "retrace, det, a and b of unit matrix" << std::endl;
  std::cout << U.retrace() << " = " << retrace(U) << std::endl << U.det()
       << " " << U.geta() << " " << U.getb() << std::endl;

  su2 U1;
  U1 = config[0] * config[1];
  std::cout << "det of U * U =" << U1.det() << std::endl;

  U1 = su2(Complex(0.8, 0.3), Complex(0.1, 0.4));
  su2 U2(0.9, Complex(0.3, 0.2));
  U1.restoreSU();
  U2.restoreSU();
  su2  U3 = U1 * U2;
  std::cout << "U3 = U1 * U2: det, retrace, a, b" << std::endl;
  std::cout << "should be: 1, 1.3264 (0.6632,0.184826) (0.293548,0.6632)" << std::endl;
  std::cout << U3.det() << " " << retrace(U1 * U2) << " " << U3.geta() << " " << U3.getb() << std::endl;

  U = config[0] * config[0].dagger();
  std::cout << "Test of U*U^dagger, should be: (1,0) (0,0)" << std::endl;
  std::cout << U.geta() << " " << U.getb() << std::endl;

  std::cout << std::endl << "Tests of U(1)" << std::endl << std::endl;;
  _u1 u;
  double a = u.retrace();
  std::cout << "test basic retrace, det and geta functions" << std::endl;
  std::cout << u.retrace() << " = " << retrace(u) << std::endl
       << u.det() << " " << u.geta() << " " << std::endl;

  _u1 x(0.5*2*M_PI), y(0.7*2*M_PI), z;
  std::cout << "test of det, the two following must be equal" << std::endl;
  std::cout << x.det() << " = " << std::exp(0.5*2*M_PI*Complex(0., 1.)) << std::endl;
  std::cout << "test of retrace, the two following must be equal" << std::endl;
  std::cout << x.retrace() << " = " << retrace(x) << std::endl;
  std::cout << "test multiplication, the two complex numbers must agree" << std::endl;

  z = x * y;
  std::cout << z.det() << " = " << x.det() * y.det() << std::endl;

  std::cout << "U(1) gauge invariance Plaquette and top. charge" << std::endl;
  
  gaugeconfig<_u1> cU(4, 4, 4, 4, 4, 1.0);

  hotstart(cU, 124665, 0.);

  double plaquette = flat_spacetime::gauge_energy(cU);
  double res = 0., Q = 0.;
  std::cout << "Initital Plaquette: " << plaquette/cU.getVolume()/6. << std::endl; 
  
  random_gauge_trafo(cU, 654321);
  plaquette = flat_spacetime::gauge_energy(cU);
  std::cout << "Plaquette after rnd trafo: " << plaquette/cU.getVolume()/6. << std::endl; 

  // set all links to 1
  hotstart(cU, 124665, 0.);
  // now a specific configuration
  // for non-zero top charge
  std::vector<size_t> xz = {1, 1, 1, 1};
  cU(xz, 0).set(0);
  cU(xz, 1).set(pi()/2.);

  cU(xz, 2).set(pi()/2);
  cU(xz, 3).set(0);

  xz = {0, 1, 1, 1};
  cU(xz, 0).set(pi());
  xz = {1, 0, 1, 1};
  cU(xz, 1).set(pi()/2.);

  xz = {1, 1, 0, 0};
  cU(xz, 2).set(pi());
  //cU(xz, 2).set(0);
  xz = {1, 1, 1, 0};
  cU(xz, 3).set(pi()/2);

  flat_spacetime::energy_density(cU, res, Q);

  std::cout << "charge: " << Q << std::endl;
  std::cout << "should be: 0.0126651" << std::endl;
  random_gauge_trafo(cU, 654321);
  Q = 0;
  flat_spacetime::energy_density(cU, res, Q);
  std::cout << "Charge after random gauge trafo: " << Q << std::endl;
  return(0);
}
