#pragma once

#include"su2.hh"
#include<vector>

using std::vector;

class gaugeconfig {
public:
  using value_type = su2;

  gaugeconfig(const size_t Ls, const size_t Lt) : 
    Ls(Ls), Lt(Lt), volume(Ls*Ls*Ls*Lt) {
    data.resize(volume*4);
  }
  size_t getLs() const {
    return(Ls);
  }
  size_t getLt() const {
    return(Lt);
  }
  size_t getVolume() const {
    return(volume);
  }
  size_t getSize() const {
    return(volume*4);
  }
  value_type &operator()(size_t const t, size_t const x, size_t const y, size_t const z, size_t const mu) {
    return data[getIndex(t, x, y, z, mu)];
  }

  const value_type &operator()(size_t const t, size_t const x, size_t const y, size_t const z, size_t const mu) const {
    return data[getIndex(t, x, y, z, mu)];
  }

  value_type &operator()(std::vector<size_t> const &coords, size_t const mu) {
    return data[getIndex(coords[0], coords[1], coords[2], coords[3], mu)];
  }

  const value_type &operator()(std::vector<size_t> const &coords, size_t const mu) const {
    return data[getIndex(coords[0], coords[1], coords[2], coords[3], mu)];
  }

  value_type &operator[](size_t const index) {
    return data[index];
  }

  const value_type &operator[](size_t const index) const {
    return data[index];
  }
  

private:
  size_t Ls, Lt, volume;

  vector<value_type> data;

  size_t getIndex(const size_t t, const size_t x, const size_t y, const size_t z, const size_t mu) const {
    size_t y0 = (t + Lt) % Lt;
    size_t y1 = (x + Ls) % Ls;
    size_t y2 = (y + Ls) % Ls;
    size_t y3 = (z + Ls) % Ls;
    size_t _mu = (mu + 4) % 4;
    return( (((y0*Ls + y1)*Ls + y2)*Ls + y3)*4 + _mu );
  }
};

gaugeconfig coldstart(size_t Ls, size_t Lt);
gaugeconfig hotstart(size_t Ls, size_t Lt, 
                     const int seed, const double delta);
