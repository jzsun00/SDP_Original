/*
  Jiazheng Sun
  Updated: Jun 19, 2024

  Class:
  FockState<NumsType>
  SpinBaseState<NumsType>
  State<StateType>
  Basis<BaseStateType>
  
  Define basis states: FockState for Fermions and Bosons,
  SpinBaseState for spin systems.
  Also define general quantum states as superposition of basis states.
  Fermi and Boson systems should inherit FockState and State.
  Spin systems should inherit SpinBaseState and State.
*/

#ifndef QM_STATES_HPP
#define QM_STATES_HPP

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>

#include "./settings.hpp"

//----------------------------------------------------------FockState<NumsType>----------

template<typename NumsType>
class FockState {
 protected:
  std::vector<NumsType> Nums;

 public:
  /*Construct a Fock state for general system.*/
  FockState() : Nums() {}
  FockState(std::vector<NumsType> & input) : Nums(input) {}
  FockState(FockState const & rhs) : Nums(rhs.Nums) {}
  ~FockState() {}
  /*Get information of the Fock state.*/
  size_t getSize() const { return Nums.size(); }
  std::vector<NumsType> getNums() const { return Nums; };
  std::string toString() const;
  /*Overload operators.*/
  FockState & operator=(FockState const & rhs);
  bool operator==(FockState const & rhs) const { return Nums == rhs.Nums; }
  bool operator!=(FockState const & rhs) const { return !(*this == rhs); }
  bool operator[](size_t n) const { return Nums[n]; }
};

//-------------------------------------------------------SpinBaseState<NumsType>---------

template<typename NumsType>
class SpinBaseState {
 protected:
  std::vector<NumsType> Nums;

 public:
  /*Construct a base state for general spin system.
    Default constructor uses an empty std::vector.*/
  SpinBaseState() : Nums() {}
  SpinBaseState(const std::vector<NumsType> & input) : Nums(input) {}
  SpinBaseState(const SpinBaseState<NumsType> & rhs) : Nums(rhs.Nums) {}
  virtual ~SpinBaseState() {}
  /*Get information of the spin base state.*/
  size_t getSize() const { return Nums.size(); }
  std::vector<NumsType> getNums() const { return Nums; };
  std::string toString() const;
  virtual std::string numToString(NumsType num) const = 0;
  /*Overload operators.*/
  SpinBaseState & operator=(const SpinBaseState & rhs);
  bool operator==(const SpinBaseState & rhs) const { return Nums == rhs.Nums; }
  bool operator[](size_t n) const { return Nums[n]; }
};

//-----------------------------------------------------------State<StateType>------------

template<typename StateType>
class State {
 protected:
  std::vector<std::pair<std::complex<double>, StateType> > Terms;

 public:
  typedef std::pair<std::complex<double>, StateType> TermType;
  /*Construct a general quantum state.
    Default constructor uses an empty std::vector.*/
  State() : Terms() {}
  State(const StateType & fs) : Terms(1, TermType(std::complex<double>(1.0, 0), fs)) {}
  State(std::complex<double> pref, const StateType & fs) : Terms(1, TermType(pref, fs)) {}
  State(const State & rhs) : Terms(rhs.Terms) {}
  virtual ~State() {}
  /*Get information of the general quantum state.*/
  size_t getSize() const { return Terms.size(); }
  typename std::vector<TermType>::const_iterator getBegin() const {
    return Terms.begin();
  }
  typename std::vector<TermType>::const_iterator getEnd() const { return Terms.end(); }
  std::string toString() const;
  /*Overload operators.*/
  std::pair<std::complex<double>, StateType> operator[](size_t n) const;
  State & operator=(const State & rhs);
  State & operator+=(const StateType & rhs);
  State & operator+=(const TermType & rhs);
  State & operator+=(const State & rhs);
  State & operator-=(const StateType & rhs);
  State & operator-=(const TermType & rhs);
  State & operator-=(const State & rhs);
  State & operator*=(std::complex<double> pref);
  void eraseZeros();

 protected:
  /*Find the same Fock state for += and -= operation.
    Return the corresponding iterator if same Fock state is found,
    otherwise return Terms.end().*/
  typename std::vector<TermType>::iterator findSameFockState(const StateType & fs);
};

//--------------------------------------------------------Basis<BaseStateType>-----------

template<typename BaseStateType>
class Basis {
 protected:
  std::vector<BaseStateType> States;

 public:
  /*Construct a basis.
    Default constructor uses an empty std::vector.*/
  Basis() : States() {}
  virtual ~Basis() {}
  /*Fill the basis.*/
  virtual void init() = 0;
  /*Get information of the basis.*/
  size_t getSize() const { return States.size(); }
  virtual std::string toString() = 0;
  typename std::vector<BaseStateType>::const_iterator getBegin() const {
    return States.begin();
  }
  typename std::vector<BaseStateType>::const_iterator getEnd() const {
    return States.end();
  }
  /*Overload operators.*/
  BaseStateType operator[](size_t n) const { return States[n]; }
};

//------------------------------------------------------------Inner Product--------------
template<typename StateType>
double innerProduct(StateType lhs, StateType rhs);

template<typename StateType>
std::complex<double> innerProduct(StateType lhs, State<StateType> rhs);

template<typename StateType>
std::complex<double> innerProduct(State<StateType> lhs, StateType rhs);

template<typename StateType>
std::complex<double> innerProduct(State<StateType> lhs, State<StateType> rhs);

#endif  //QM_STATES_HPP
