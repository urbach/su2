#pragma once
#include"hamiltonian_field.hh"
#include<vector>

// mother class for monomials
template<class T> class monomial {
public:
  double Hold, Hnew;
  monomial(unsigned int _timescale) : Hold(0.), Hnew(0.), timescale(_timescale), mdactive(true) {}
  monomial() : Hold(0.), Hnew(0.), timescale(0), mdactive(true) {}
  virtual void heatbath(hamiltonian_field<T> const &h) = 0;
  virtual void accept(hamiltonian_field<T> const &h) = 0;
  virtual void derivative(adjointfield<T> &deriv, hamiltonian_field<T> const &h, const T fac) const = 0;
  void reset() {
    Hold = 0.;
    Hnew = 0.;
  }
  void setTimescale(unsigned int _timescale) {
    timescale = _timescale;
  }
  unsigned int getTimescale() const {
    return timescale;
  }
  double getDeltaH() const {
    return Hnew - Hold;
  }
  void setmdactive() {
    mdactive = true;
  }
  void setmdpassive() {
    mdactive = false;
  }
  bool getmdactive() const {
    return mdactive;
  }
private:
  unsigned int timescale;
  bool mdactive;
};

template<class T> class kineticmonomial : public monomial<T> {
public:
  kineticmonomial(unsigned int _timescale) : monomial<T>::monomial(_timescale) {}
  void heatbath(hamiltonian_field<T> const &h) override {
    monomial<T>::Hold = 0.5 * ((*h.momenta)*(*h.momenta));
  }
  void accept(hamiltonian_field<T> const &h) override {
    monomial<T>::Hnew = 0.5 * ((*h.momenta)*(*h.momenta));
  }
  virtual void derivative(adjointfield<T> &deriv, hamiltonian_field<T> const &h, const T fac) const {}
};

#include"gaugemonomial.hh"
