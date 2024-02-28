/*
  Jiazheng Sun
  Updated: Feb 27, 2024
  Define Fock states, general quantum states for Fermionic systems.
  Define the full basis for a Fermion lattice system.
*/

#ifndef ORI_SDP_GS_FERMISTATES_HPP
#define ORI_SDP_GS_FERMISTATES_HPP

#include <bitset>
#include <cmath>
#include <numeric>

#include "operators.hpp"

using std::vector;

//----------------------------------------------------------------FermiFockState--------

class FermiFockState {
 protected:
  vector<bool> Nums;

 public:
  /*Construct a Fock state for Fermion system.*/
  FermiFockState() : Nums() {}
  FermiFockState(size_t n) : Nums(n, false) {}
  FermiFockState(vector<bool> & input) : Nums(input) {}
  FermiFockState(FermiFockState const & rhs) : Nums(rhs.Nums) {}
  ~FermiFockState() {}
  /*Get information of the Fock state.*/
  size_t getSize() const { return Nums.size(); }
  vector<bool> getNums() const { return Nums; };
  std::string toString();
  /*Overload operators.*/
  FermiFockState & operator=(FermiFockState const & rhs);
  bool operator==(FermiFockState const & rhs) const { return Nums == rhs.Nums; }
  bool operator[](size_t n) const { return Nums[n]; }
};

//-------------------------------------------------------------------FermiState---------

class FermiState {
 protected:
  vector<pair<complex<double>, FermiFockState> > Terms;

 public:
  typedef pair<complex<double>, FermiFockState> TermType;
  /*Construct a general quantum state for Fermions.*/
  FermiState() : Terms() {}
  FermiState(FermiFockState const & ffs) : Terms(1) {
    Terms[0].first = complex<double>(1, 0);
    Terms[0].second = ffs;
  }
  FermiState(complex<double> pref, FermiFockState const & ffs) : Terms(1) {
    Terms[0].first = pref;
    Terms[0].second = ffs;
  }
  FermiState(FermiState const & rhs) : Terms(rhs.Terms) {}
  ~FermiState() {}
  /*Get information of the general quantum state.*/
  size_t getSize() const { return Terms.size(); }
  vector<pair<complex<double>, FermiFockState> >::iterator getBegin() {
    return Terms.begin();
  }
  vector<pair<complex<double>, FermiFockState> >::iterator getEnd() {
    return Terms.end();
  }
  std::string toString();
  /*Overload operators.*/
  pair<complex<double>, FermiFockState> operator[](size_t n) const;
  FermiState & operator+=(FermiFockState const & rhs);
  FermiState & operator+=(TermType const & rhs);
  FermiState & operator+=(FermiState & rhs);
  FermiState & operator-=(FermiFockState const & rhs);
  FermiState & operator-=(TermType const & rhs);
  FermiState & operator-=(FermiState & rhs);
  void eraseZeros();

 private:
  /*Find the same Fock state for += operation.
    Return the corresponding iterator if same Fock state is found,
    otherwise return Terms.end().*/
  vector<pair<complex<double>, FermiFockState> >::iterator findSameFockState(
      FermiFockState const & ffs);
};

//-------------------------------------------------------------------FermiBasis---------

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
  /*Overload operators.*/
  FermiFockState operator[](size_t n) const { return States[n]; }
};

//----------------------------------------------------------------Other Functions-------

double innerProduct(FermiFockState lhs, FermiFockState rhs);
complex<double> innerProduct(FermiFockState lhs, FermiState rhs);
complex<double> innerProduct(FermiState lhs, FermiFockState rhs);
complex<double> innerProduct(FermiState lhs, FermiState rhs);

#endif
