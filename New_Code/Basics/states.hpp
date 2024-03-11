/*
  Jiazheng Sun
  Updated: Mar 10, 2024

  Define Fock states and quantum states for general systems.
*/

#ifndef ORI_SDP_GS_STATES_HPP
#define ORI_SDP_GS_STATES_HPP

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

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

//-----------------------------------------------------------------------State-----------

template<typename StateType>
class State {
 protected:
  vector<pair<complex<double>, StateType> > Terms;

 public:
  typedef pair<complex<double>, StateType> TermType;
  /*Construct a general quantum state for Fermions.*/
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

//-------------------------------------------------------------------Inner Product-------
template<typename StateType>
double innerProduct(StateType lhs, StateType rhs);

template<typename StateType>
complex<double> innerProduct(StateType lhs, State<StateType> rhs);

template<typename StateType>
complex<double> innerProduct(State<StateType> lhs, StateType rhs);

template<typename StateType>
complex<double> innerProduct(State<StateType> lhs, State<StateType> rhs);

#include "states.cpp"

#endif  //ORI_SDP_GS_STATES_HPP
