/*
  Jan 31, 2024  Jiazheng Sun

  Define ladder operators, monomials and polynomials.
*/

#ifndef ORI_SDP_GS_FERMISTATES_HPP
#define ORI_SDP_GS_FERMISTATES_HPP

#include <bitset>
#include <cmath>

#include "operators.hpp"

using std::vector;

class FermiFockState {
 protected:
  vector<bool> Nums;

 public:
  /*Construct a Fock state for Fermion system.*/
  FermiFockState() : Nums() {}
  FermiFockState(size_t n) : Nums(n, false) {}
  FermiFockState(vector<bool> & input) : Nums(input) {}
  //FermiFockState(FermiFockState & rhs) : States(rhs.States) {}
  ~FermiFockState() {}
  size_t getSize() const { return Nums.size(); }
  std::string toString();
  bool operator==(FermiFockState const & rhs) const;
  bool operator[](size_t n) const { return Nums[n]; }
};

class FermiState {
 protected:
  vector<pair<complex<double>, FermiFockState> > Terms;

 public:
  typedef pair<complex<double>, FermiFockState> TermType;
  FermiState() : Terms() {}
  FermiState(FermiFockState const & ffs) : Terms(1) {
    Terms[0].first = complex<double>(1, 0);
    Terms[0].second = ffs;
  }
  size_t getSize() const { return Terms.size(); }
  vector<pair<complex<double>, FermiFockState> >::iterator getBegin() {
    return Terms.begin();
  }
  vector<pair<complex<double>, FermiFockState> >::iterator getEnd() {
    return Terms.end();
  }
  std::string toString();
  pair<complex<double>, FermiFockState> operator[](size_t n) const { return Terms[n]; }
  FermiState & operator+=(TermType const & rhs);
  FermiState & operator+=(FermiState & rhs);
  FermiState & operator-=(TermType const & rhs);
  FermiState & operator-=(FermiState & rhs);
};

class FermiBasis {
 protected:
  size_t Sites;
  vector<FermiFockState> States;

 public:
  FermiBasis() : Sites(0), States() {}
  FermiBasis(size_t n) : Sites(n), States() {}
  void init();
  std::string toString();
  vector<FermiFockState>::iterator getBegin() { return States.begin(); }
  vector<FermiFockState>::iterator getEnd() { return States.end(); }
  FermiFockState operator[](size_t n) const { return States[n]; }
};

double innerProduct(FermiFockState lhs, FermiFockState rhs);
complex<double> innerProduct(FermiFockState lhs, FermiState rhs);
complex<double> innerProduct(FermiState lhs, FermiFockState rhs);
complex<double> innerProduct(FermiState lhs, FermiState rhs);

#endif
