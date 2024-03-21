/*
  Jiazheng Sun
  Updated: Mar 19, 2024

  Define Fock states, quantum states for Fermionic systems.
  Define the full basis for a Fermion lattice system.
*/

#ifndef ORI_SDP_GS_FERMISTATES_HPP
#define ORI_SDP_GS_FERMISTATES_HPP

#include <bitset>
#include <cmath>
#include <numeric>

#include "../Basics/states.hpp"

//----------------------------------------------------------------FermiFockState--------

class FermiFockState : public FockState<bool> {
 public:
  /*Construct a Fock state for Fermion system.
    Constructors are identical to FockState.*/
  FermiFockState() : FockState() {}
  FermiFockState(vector<bool> & input) : FockState(input) {}
  FermiFockState(FockState<bool> const & rhs) : FockState(rhs) {}
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
  FermiState(FermiState const & rhs) : State(rhs) {}
  ~FermiState() {}
};

//-------------------------------------------------------------------FermiBasis---------

class Fermi1DBasis : Basis<FermiFockState> {
 protected:
  size_t Sites;

 public:
  /*Construct the entire basis of Fermi systems.*/
  Fermi1DBasis() : Basis<FermiFockState>(), Sites(0) {}
  Fermi1DBasis(size_t n) : Basis<FermiFockState>(), Sites(n) {}
  void init();
  /*Get information of the full basis.*/
  std::string toString();
};

#endif
