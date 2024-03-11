/*
  Jiazheng Sun
  Updated: Mar 10, 2024

  Define Fock states, quantum states for spin systems.
  Define the full basis for a spin lattice system.
*/

#ifndef ORI_SDP_GS_SPINSTATES_HPP
#define ORI_SDP_GS_SPINSTATES_HPP

#include <bitset>
#include <cmath>
#include <numeric>

#include "operators.hpp"
#include "states.hpp"

//----------------------------------------------------------------SpinBaseState----------

template<typename NumsType>
class SpinBaseState {
 protected:
  vector<NumsType> Nums;

 public:
  /*Construct a pure state for general spin system.*/
  SpinBaseState() : Nums() {}
  SpinBaseState(vector<NumsType> & input) : Nums(input) {}
  SpinBaseState(SpinBaseState const & rhs) : Nums(rhs.Nums) {}
  ~SpinBaseState() {}
  /*Get information of the Fock state.*/
  size_t getSize() const { return Nums.size(); }
  vector<NumsType> getNums() const { return Nums; };
  std::string toString() const;
  /*Overload operators.*/
  SpinBaseState & operator=(SpinBaseState const & rhs);
  bool operator==(SpinBaseState const & rhs) const { return Nums == rhs.Nums; }
  bool operator[](size_t n) const { return Nums[n]; }
};

//-------------------------------------------------------------SpinHalfBaseState--------

class SpinHalfBaseState : public SpinBaseState<bool> {
 public:
  SpinHalfBaseState() : SpinBaseState() {}
  SpinHalfBaseState(vector<bool> & input) : SpinBaseState(input) {}
  SpinHalfBaseState(SpinBaseState const & rhs) : SpinBaseState(rhs) {}
  ~SpinHalfBaseState() {}
};

//------------------------------------------------------------------SpinState-----------

template<typename SpinType>
class SpinState : public State<SpinType> {
 public:
  typedef pair<complex<double>, SpinHalfBaseState> TermType;
  /*Construct a general quantum state for Fermions.
    Constructors are identical to State*/
  SpinState() : State<SpinType>() {}
  SpinState(SpinHalfBaseState const & ffs) : State<SpinType>(ffs) {}
  SpinState(complex<double> pref, SpinHalfBaseState const & ffs) :
      State<SpinType>(pref, ffs) {}
  SpinState(State<SpinType> const & rhs) : State<SpinType>(rhs) {}
  ~SpinState() {}
};

//----------------------------------------------------------------SpinHalfState---------

class SpinHalfState : public SpinState<SpinHalfBaseState> {
 public:
  typedef pair<complex<double>, SpinHalfBaseState> TermType;
  /*Construct a general quantum state for Fermions.
    Constructors are identical to State*/
  SpinHalfState() : SpinState() {}
  SpinHalfState(SpinHalfBaseState const & ffs) : SpinState(ffs) {}
  SpinHalfState(complex<double> pref, SpinHalfBaseState const & ffs) :
      SpinState(pref, ffs) {}
  SpinHalfState(State const & rhs) : SpinState(rhs) {}
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

#include "spinStates.cpp"

#endif
