/*
  Jiazheng Sun
  Updated: Mar 12, 2024

  Define Fock states, quantum states for spin systems.
  Define the full basis for a spin lattice system.
*/

#ifndef ORI_SDP_GS_SPINSTATES1D_HPP
#define ORI_SDP_GS_SPINSTATES1D_HPP

#include <bitset>
#include <cmath>
#include <numeric>

#include "../Basics/states.hpp"

//-------------------------------------------------------------SpinHalfBaseState--------

class SpinHalfBaseState : public SpinBaseState<bool> {
 public:
  SpinHalfBaseState() : SpinBaseState() {}
  SpinHalfBaseState(size_t N) : SpinBaseState() { Nums = vector<bool>(N, false); }
  SpinHalfBaseState(vector<bool> & input) : SpinBaseState(input) {}
  SpinHalfBaseState(SpinBaseState const & rhs) : SpinBaseState(rhs) {}
  ~SpinHalfBaseState() {}
};

//----------------------------------------------------------------SpinHalfState---------

class SpinHalfState : public State<SpinHalfBaseState> {
 public:
  typedef pair<complex<double>, SpinHalfBaseState> TermType;
  /*Construct a general quantum state for Fermions.
    Constructors are identical to State*/
  SpinHalfState() : State() {}
  SpinHalfState(SpinHalfBaseState const & ffs) : State(ffs) {}
  SpinHalfState(complex<double> pref, SpinHalfBaseState const & ffs) : State(pref, ffs) {}
  SpinHalfState(State const & rhs) : State(rhs) {}
  ~SpinHalfState() {}
};

//-------------------------------------------------------------------SpinHalfBasis-------

class SpinHalfBasis {
 protected:
  size_t Sites;
  vector<SpinHalfBaseState> States;

 public:
  /*Construct the entire basis of Fermi systems.*/
  SpinHalfBasis() : Sites(0), States() {}
  SpinHalfBasis(size_t n) : Sites(n), States() {}
  void init();
  /*Get information of the full basis.*/
  std::string toString();
  vector<SpinHalfBaseState>::const_iterator getBegin() const { return States.begin(); }
  vector<SpinHalfBaseState>::const_iterator getEnd() const { return States.end(); }
  /*Overload operators.*/
  SpinHalfBaseState operator[](size_t n) const { return States[n]; }
};

//#include "spinStates1D_Tem.cpp"

#endif
