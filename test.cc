#include"su2.hh"
#include"u1.hh"
#include"random_element.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"parse_commandline.hh"
#include"energy_density.hh"

#include<iostream>
#include<vector>

using std::vector;
using std::cout;
using std::endl;


int main() {
  size_t Lx = 8, Ly = 8, Lz = 8, Lt = 16;


  cout << "Tests of SU(2)" << endl << endl;;
  
  vector<su2> config;
  config.resize(Lx*Ly*Lz*Lt);
  
  // set all to 1
  for(auto i = config.begin(), end =config.end(); i < end; i++) {
    *i = su2(1., 0.);
  }

  su2 U(1., 0.);
  cout << "retrace, det, a and b of unit matrix" << endl;
  cout << U.retrace() << " = " << retrace(U) << endl << U.det()
       << " " << U.geta() << " " << U.getb() << endl;

  su2 U1;
  U1 = config[0] * config[1];
  cout << "det of U * U =" << U1.det() << endl;

  U1 = su2(Complex(0.8, 0.3), Complex(0.1, 0.4));
  su2 U2(0.9, Complex(0.3, 0.2));
  U1.restoreSU();
  U2.restoreSU();
  su2  U3 = U1 * U2;
  cout << "U3 = U1 * U2: det, retrace, a, b" << endl;
  cout << "should be: 1, 1.3264 (0.6632,0.184826) (0.293548,0.6632)" << endl;
  cout << U3.det() << " " << retrace(U1 * U2) << " " << U3.geta() << " " << U3.getb() << endl;

  U = config[0] * config[0].dagger();
  cout << "Test of U*U^dagger, should be: (1,0) (0,0)" << endl;
  cout << U.geta() << " " << U.getb() << endl;

  cout << endl << "Tests of U(1)" << endl << endl;;
  _u1 u;
  double a = u.retrace();
  cout << "test basic retrace, det and geta functions" << endl;
  cout << u.retrace() << " = " << retrace(u) << endl
       << u.det() << " " << u.geta() << " " << endl;

  _u1 x(0.5*2*M_PI), y(0.7*2*M_PI), z;
  cout << "test of det, the two following must be equal" << endl;
  cout << x.det() << " = " << std::exp(0.5*2*M_PI*Complex(0., 1.)) << endl;
  cout << "test of retrace, the two following must be equal" << endl;
  cout << x.retrace() << " = " << retrace(x) << endl;
  cout << "test multiplication, the two complex numbers must agree" << endl;

  z = x * y;
  cout << z.det() << " = " << x.det() * y.det() << endl;

  cout << "U(1) gauge invariance Plaquette and top. charge" << endl;
  
  gaugeconfig<_u1> cU(4, 4, 4, 4, 4, 1.0);

  hotstart(cU, 124665, 0.);

  double plaquette = gauge_energy(cU);
  double res = 0., Q = 0.;
  cout << "Initital Plaquette: " << plaquette/cU.getVolume()/6. << endl; 
  
  random_gauge_trafo(cU, 654321);
  plaquette = gauge_energy(cU);
  cout << "Plaquette after rnd trafo: " << plaquette/cU.getVolume()/6. << endl; 

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

  energy_density(cU, res, Q);

  cout << "charge: " << Q << endl;
  cout << "should be: 0.0126651" << endl;
  random_gauge_trafo(cU, 654321);
  Q = 0;
  energy_density(cU, res, Q);
  cout << "Charge after random gauge trafo: " << Q << endl;
  return(0);
}
