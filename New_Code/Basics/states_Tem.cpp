/*
  Jiazheng Sun
  Updated: Mar 12, 2024

  Implementations of methods in class:
  Fockstate, State.
  Also implementations of inner product.
 */

#ifndef ORI_SDP_GS_STATES_TEM_CPP
#define ORI_SDP_GS_STATES_TEM_CPP

#include "states.hpp"

//-------------------------------------------------------------------FockState-----------

template<typename NumsType>
std::string FockState<NumsType>::toString() const {
  std::string ans = "|";
  for (typename vector<NumsType>::const_iterator it = Nums.begin(); it != Nums.end();
       ++it) {
    ans += " ";
    ans += std::to_string(*it);
    ans += ",";
  }
  ans.pop_back();
  ans += " >";
  return ans;
}

template<typename NumsType>
FockState<NumsType> & FockState<NumsType>::operator=(FockState<NumsType> const & rhs) {
  Nums = rhs.Nums;
  return *this;
}

//----------------------------------------------------------------SpinBaseState----------

template<typename NumsType>
std::string SpinBaseState<NumsType>::toString() const {
  std::string ans = "|";
  for (typename vector<NumsType>::const_iterator it = Nums.begin(); it != Nums.end();
       ++it) {
    ans += " ";
    ans += std::to_string(*it);
    ans += ",";
  }
  ans.pop_back();
  ans += " >";
  return ans;
}

template<typename NumsType>
SpinBaseState<NumsType> & SpinBaseState<NumsType>::operator=(
    SpinBaseState<NumsType> const & rhs) {
  Nums = rhs.Nums;
  return *this;
}

//-----------------------------------------------------------------------State-----------

template<typename StateType>
std::string State<StateType>::toString() const {
  std::string ans = "";
  for (typename vector<pair<complex<double>, StateType> >::const_iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    ans += "  (";
    ans += std::to_string(it->first.real());
    ans += " + ";
    ans += std::to_string(it->first.imag());
    ans += ")";
    ans += it->second.toString();
    ans += "  +\n";
  }
  ans.pop_back();
  ans.pop_back();
  return ans;
}

template<typename StateType>
pair<complex<double>, StateType> State<StateType>::operator[](size_t n) const {
  return Terms.at(n);
}

template<typename StateType>
State<StateType> & State<StateType>::operator=(State<StateType> const & rhs) {
  Terms = rhs.Terms;
  return *this;
}

template<typename StateType>
typename vector<pair<complex<double>, StateType> >::iterator
State<StateType>::findSameFockState(StateType const & ffs) {
  for (typename vector<pair<complex<double>, StateType> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    if (it->second == ffs) {
      return it;
    }
  }
  return Terms.end();
}

//+=
template<typename StateType>
State<StateType> & State<StateType>::operator+=(StateType const & rhs) {
  pair<complex<double>, StateType> toAdd(complex<double>(1, 0), rhs);
  *this += toAdd;
  return *this;
}

template<typename StateType>
State<StateType> & State<StateType>::operator+=(
    pair<complex<double>, StateType> const & rhs) {
  if (std::abs(rhs.first) < ERROR) {
    return *this;
  }
  typename vector<pair<complex<double>, StateType> >::iterator it =
      findSameFockState(rhs.second);
  if (it == Terms.end()) {
    Terms.push_back(rhs);
  }
  else {
    it->first += rhs.first;
  }
  return *this;
}

template<typename StateType>
State<StateType> & State<StateType>::operator+=(State<StateType> const & rhs) {
  for (typename vector<pair<complex<double>, StateType> >::const_iterator termIt =
           rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this += *termIt;
  }
  return *this;
}

//-=
template<typename StateType>
State<StateType> & State<StateType>::operator-=(StateType const & rhs) {
  pair<complex<double>, StateType> toAdd(complex<double>(-1, 0), rhs);
  *this += toAdd;
  return *this;
}

template<typename StateType>
State<StateType> & State<StateType>::operator-=(
    pair<complex<double>, StateType> const & rhs) {
  if (std::abs(rhs.first) < ERROR) {
    return *this;
  }
  typename vector<pair<complex<double>, StateType> >::iterator it =
      findSameFockState(rhs.second);
  if (it == Terms.end()) {
    pair<complex<double>, StateType> copy(complex<double>(-1, 0) * rhs.first, rhs.second);
    Terms.push_back(copy);
  }
  else {
    it->first -= rhs.first;
  }
  return *this;
}

template<typename StateType>
State<StateType> & State<StateType>::operator-=(State<StateType> const & rhs) {
  for (typename vector<pair<complex<double>, StateType> >::const_iterator termIt =
           rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this -= *termIt;
  }
  return *this;
}

//*=
template<typename StateType>
State<StateType> & State<StateType>::operator*=(complex<double> pref) {
  for (typename vector<pair<complex<double>, StateType> >::iterator termIt =
           Terms.begin();
       termIt != Terms.end();
       ++termIt) {
    termIt->first *= pref;
  }
  return *this;
}

template<typename StateType>
bool isZeroS(pair<complex<double>, StateType> term) {
  return std::abs(term.first) < ERROR;
}

template<typename StateType>
void State<StateType>::eraseZeros() {
  Terms.erase(std::remove_if(Terms.begin(), Terms.end(), isZeroS), Terms.end());
}

//-------------------------------------------------------------------Inner Product-------

template<typename StateType>
double innerProduct(StateType lhs, StateType rhs) {
  if (lhs == rhs) {
    return 1;
  }
  return 0;
}

template<typename StateType>
complex<double> innerProduct(StateType lhs, State<StateType> rhs) {
  complex<double> ans(0, 0);
  for (typename vector<pair<complex<double>, StateType> >::const_iterator termIt =
           rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    ans += termIt->first * innerProduct(lhs, termIt->second);
  }
  return ans;
}

template<typename StateType>
complex<double> innerProduct(State<StateType> lhs, StateType rhs) {
  return std::conj(innerProduct(rhs, lhs));
}

template<typename StateType>
complex<double> innerProduct(State<StateType> lhs, State<StateType> rhs) {
  complex<double> ans(0, 0);
  for (typename vector<pair<complex<double>, StateType> >::const_iterator it =
           lhs.getBegin();
       it != lhs.getEnd();
       ++it) {
    ans += it->first * innerProduct(it->second, rhs);
  }
  return ans;
}

#endif  //ORI_SDP_GS_STATES_TEM_CPP
