/*
  Jiazheng Sun
  Updated: Mar 18, 2024

  Class:
  FockState, SpinBaseState, State, Basis.

  Define basis states: FockState for Fermions and Bosons,
  SpinBaseState for spin systems.
  Also define general quantum states.
  Fermi and Boson systems should inherit FockState and State.
  Spin systems should inherit SpinBaseState and State.

  Implementations for all classes are in states_Tem.cpp.
*/

#ifndef ORI_SDP_GS_STATES_HPP
#define ORI_SDP_GS_STATES_HPP

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "./settings.hpp"

using std::complex;
using std::pair;
using std::vector;

//-------------------------------------------------------------------FockState-----------

template<typename NumsType>
class FockState {
 protected:
  vector<NumsType> Nums;

 public:
  /*Construct a Fock state for general system.*/
  FockState() : Nums() {}
  FockState(vector<NumsType> & input) : Nums(input) {}
  FockState(FockState const & rhs) : Nums(rhs.Nums) {}
  ~FockState() {}
  /*Get information of the Fock state.*/
  size_t getSize() const { return Nums.size(); }
  vector<NumsType> getNums() const { return Nums; };
  std::string toString() const;
  /*Overload operators.*/
  FockState & operator=(FockState const & rhs);
  bool operator==(FockState const & rhs) const { return Nums == rhs.Nums; }
  bool operator[](size_t n) const { return Nums[n]; }
};

//----------------------------------------------------------------SpinBaseState----------

template<typename NumsType>
class SpinBaseState {
 protected:
  vector<NumsType> Nums;

 public:
  /*Construct a pure base state for general spin system.*/
  SpinBaseState() : Nums() {}
  SpinBaseState(vector<NumsType> & input) : Nums(input) {}
  SpinBaseState(SpinBaseState const & rhs) : Nums(rhs.Nums) {}
  ~SpinBaseState() {}
  /*Get information of the spin base state.*/
  size_t getSize() const { return Nums.size(); }
  vector<NumsType> getNums() const { return Nums; };
  std::string toString() const;
  /*Overload operators.*/
  SpinBaseState & operator=(SpinBaseState const & rhs);
  bool operator==(SpinBaseState const & rhs) const { return Nums == rhs.Nums; }
  bool operator[](size_t n) const { return Nums[n]; }
};

//-----------------------------------------------------------------------State-----------

template<typename StateType>
class State {
 protected:
  vector<pair<complex<double>, StateType> > Terms;

 public:
  typedef pair<complex<double>, StateType> TermType;
  /*Construct a general quantum state.*/
  State() : Terms() {}
  State(StateType const & fs) : Terms(1) {
    Terms[0].first = complex<double>(1, 0);
    Terms[0].second = fs;
  }
  State(complex<double> pref, StateType const & fs) : Terms(1) {
    Terms[0].first = pref;
    Terms[0].second = fs;
  }
  State(State const & rhs) : Terms(rhs.Terms) {}
  ~State() {}
  /*Get information of the general quantum state.*/
  size_t getSize() const { return Terms.size(); }
  typename vector<pair<complex<double>, StateType> >::const_iterator getBegin() const {
    return Terms.begin();
  }
  typename vector<pair<complex<double>, StateType> >::const_iterator getEnd() const {
    return Terms.end();
  }
  std::string toString() const;
  /*Overload operators.*/
  pair<complex<double>, StateType> operator[](size_t n) const;
  State & operator=(State const & rhs);
  State & operator+=(StateType const & rhs);
  State & operator+=(TermType const & rhs);
  State & operator+=(State const & rhs);
  State & operator-=(StateType const & rhs);
  State & operator-=(TermType const & rhs);
  State & operator-=(State const & rhs);
  State & operator*=(complex<double> pref);
  void eraseZeros();

 protected:
  /*Find the same Fock state for += operation.
    Return the corresponding iterator if same Fock state is found,
    otherwise return Terms.end().*/
  typename vector<pair<complex<double>, StateType> >::iterator findSameFockState(
      StateType const & fs);
};

//----------------------------------------------------------------------Basis------------

template<typename BaseStateType>
class Basis {
 protected:
  vector<BaseStateType> States;

 public:
  /*Construct a basis, default use empty vector.*/
  Basis() : States() {}
  ~Basis() {}
  virtual void init() = 0;
  /*Get information of the basis.*/
  virtual std::string toString() = 0;
  typename vector<BaseStateType>::const_iterator getBegin() const {
    return States.begin();
  }
  typename vector<BaseStateType>::const_iterator getEnd() const { return States.end(); }
  BaseStateType operator[](size_t n) const { return States[n]; }
};

//-------------------------------------------------------------------Inner Product-------
template<typename StateType>
double innerProduct(StateType lhs, StateType rhs);

template<typename StateType>
complex<double> innerProduct(StateType lhs, State<StateType> rhs);

template<typename StateType>
complex<double> innerProduct(State<StateType> lhs, StateType rhs);

template<typename StateType>
complex<double> innerProduct(State<StateType> lhs, State<StateType> rhs);

#include "states_Tem.cpp"

#endif  //ORI_SDP_GS_STATES_HPP
