/*
  Jan 31, 2024  Jiazheng Sun

  Define ladder operators, monomials and polynomials.
*/

#ifndef ORI_SDP_GS_FERMISTATES_HPP
#define ORI_SDP_GS_FERMISTATES_HPP

#include "operators.hpp"

#include <cmath>

using std::vector;

class FermiFockState {
 protected:
  vector<bool> States;

 public:
  FermiFockState() : States() {}
  FermiFockState(size_t n) : States(n, false) {}
  FermiFockState(vector<bool> & input) : States(input) {}
  //FermiFockState(FermiFockState & rhs) : States(rhs.States) {}
  ~FermiFockState() {}
  size_t getSize() const { return States.size(); }
  std::string toString();
  bool operator==(FermiFockState const & rhs);
  bool operator[](size_t n) const { return States[n]; }
};


class FermiState {
 protected:
  vector<pair<complex<double>, FermiFockState> > Terms;

 public:
  size_t getSize() const { return Terms.size(); }
  vector<pair<complex<double>, FermiFockState> >::iterator getBegin() { return Terms.begin(); }
  vector<pair<complex<double>, FermiFockState> >::iterator getEnd() { return Terms.end(); }
  pair<complex<double>, FermiFockState> operator[](size_t n) const { return Terms[n]; }
};


class FermiBasis {
protected:
  vector<FermiFockState> States;

public:
  FermiBasis() : States() {}
  FermiBasis(size_t n) : States(std::pow(2, n)) {
    
  }
};

double innerProduct(FermiFockState lhs, FermiFockState rhs);
complex<double> innerProduct(FermiFockState lhs, FermiState rhs);
complex<double> innerProduct(FermiState lhs, FermiFockState rhs);
complex<double> innerProduct(FermiState lhs, FermiState rhs);




#endif
