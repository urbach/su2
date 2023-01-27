#pragma once
#include "hamiltonian_field.hh"
#include <vector>

// mother class for monomials
template <typename Float, class Group> class monomial {
public:
  double Hold, Hnew;
  monomial(unsigned int _timescale)
    : Hold(0.), Hnew(0.), timescale(_timescale), mdactive(true) {}
  monomial() : Hold(0.), Hnew(0.), timescale(0), mdactive(true) {}
  virtual void heatbath(hamiltonian_field<Float, Group> const &h) = 0;
  virtual void accept(hamiltonian_field<Float, Group> const &h) = 0;
  virtual void derivative(adjointfield<Float, Group> &deriv,
                          hamiltonian_field<Float, Group> const &h,
                          const Float fac) const = 0;
  void reset() {
    Hold = 0.;
    Hnew = 0.;
  }
  void setTimescale(unsigned int _timescale) { timescale = _timescale; }
  unsigned int getTimescale() const { return timescale; }
  double getDeltaH() const { return Hnew - Hold; }
  void setmdactive() { mdactive = true; }
  void setmdpassive() { mdactive = false; }
  bool getmdactive() const { return mdactive; }

private:
  unsigned int timescale;
  bool mdactive;
};

template <typename Float, class Group>
class kineticmonomial : public monomial<Float, Group> {
public:
  kineticmonomial(unsigned int _timescale)
    : monomial<Float, Group>::monomial(_timescale) {}
  void heatbath(hamiltonian_field<Float, Group> const &h) override {
    monomial<Float, Group>::Hold = 0.5 * ((*h.momenta) * (*h.momenta));
  }
  void accept(hamiltonian_field<Float, Group> const &h) override {
    monomial<Float, Group>::Hnew = 0.5 * ((*h.momenta) * (*h.momenta));
  }
  virtual void derivative(adjointfield<Float, Group> &deriv,
                          hamiltonian_field<Float, Group> const &h,
                          const Float fac) const {}
};

#include "gaugemonomial.hh"
