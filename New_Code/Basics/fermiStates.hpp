/*
  Jiazheng Sun
  Updated: Mar 10, 2024

  Define Fock states, quantum states for Fermionic systems.
  Define the full basis for a Fermion lattice system.
*/

#ifndef ORI_SDP_GS_FERMISTATES_HPP
#define ORI_SDP_GS_FERMISTATES_HPP

#include <bitset>
#include <cmath>
#include <numeric>

#include "operators.hpp"
#include "states.hpp"

//----------------------------------------------------------------FermiFockState--------

class FermiFockState : public FockState<bool> {
 public:
  /*Construct a Fock state for Fermion system.
    Constructors are identical to FockState.*/
  FermiFockState() : FockState() {}
  FermiFockState(vector<bool> & input) : FockState(input) {}
  FermiFockState(FockState const & rhs) : FockState(rhs) {}
  ~FermiFockState() {}
};

//-------------------------------------------------------------------FermiState---------

class FermiState : public State<FermiFockState> {
 public:
  typedef pair<complex<double>, FermiFockState> TermType;
  /*Construct a general quantum state for Fermions.
    Constructors are identical to State*/
  FermiState() : State() {}
  FermiState(FermiFockState const & ffs) : State(ffs) {}
  FermiState(complex<double> pref, FermiFockState const & ffs) : State(pref, ffs) {}
  FermiState(State const & rhs) : State(rhs) {}
  ~FermiState() {}
};

//-------------------------------------------------------------------FermiBasis---------

class FermiBasis {
 protected:
  size_t Sites;
  vector<FermiFockState> States;

 public:
  /*Construct the entire basis of Fermi systems.*/
  FermiBasis() : Sites(0), States() {}
  FermiBasis(size_t n) : Sites(n), States() {}
  void init();
  /*Get information of the full basis.*/
  std::string toString();
  vector<FermiFockState>::const_iterator getBegin() const { return States.begin(); }
  vector<FermiFockState>::const_iterator getEnd() const { return States.end(); }
  /*Overload operators.*/
  FermiFockState operator[](size_t n) const { return States[n]; }
};

#endif
