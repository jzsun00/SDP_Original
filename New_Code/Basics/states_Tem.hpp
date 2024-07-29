/*
  Jiazheng Sun
  Updated: Jul 28, 2024
  
  Class Implementations:
  FockState<NumsType>
  SpinBaseState<NumsType>
  State<StateType>
  
  Function Implementations:
  double innerProduct(StateType lhs, StateType rhs);
  complex<double> innerProduct(StateType lhs, State<StateType> rhs);
  complex<double> innerProduct(State<StateType> lhs, StateType rhs);
  complex<double> innerProduct(State<StateType> lhs, State<StateType> rhs);
*/

#ifndef QM_STATES_TEM_HPP
#define QM_STATES_TEM_HPP

#include "settings.hpp"
#include "states.hpp"

//----------------------------------------------------------FockState<NumsType>----------

template<typename NumsType>
std::string FockState<NumsType>::toString() const {
  std::string ans = "|";
  for (typename std::vector<NumsType>::const_iterator it = Nums.begin(); it != Nums.end();
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

//-------------------------------------------------------SpinBaseState<NumsType>---------

template<typename NumsType>
std::string SpinBaseState<NumsType>::toString() const {
  std::string ans = "|";
  for (auto it = Nums.begin(); it != Nums.end(); ++it) {
    ans += " ";
    ans += numToString(*it);
    ans += ",";
  }
  ans.pop_back();
  ans += " >";
  return ans;
}

template<typename NumsType>
SpinBaseState<NumsType> & SpinBaseState<NumsType>::operator=(
    const SpinBaseState<NumsType> & rhs) {
  Nums = rhs.Nums;
  return *this;
}

//-----------------------------------------------------------State<StateType>------------

template<typename StateType>
std::string State<StateType>::toString() const {
  std::string ans = "";
  for (auto it = Terms.begin(); it != Terms.end(); ++it) {
    ans += "  (";
    ans += complex_toString(it->first);
    ans += ") ";
    ans += it->second.toString();
    ans += "  +\n";
  }
  ans.pop_back();
  ans.pop_back();
  return ans;
}

template<typename StateType>
std::pair<std::complex<double>, StateType> State<StateType>::operator[](size_t n) const {
  return Terms[n];
}

template<typename StateType>
State<StateType> & State<StateType>::operator=(State<StateType> const & rhs) {
  Terms = rhs.Terms;
  return *this;
}

template<typename StateType>
typename std::vector<std::pair<std::complex<double>, StateType> >::iterator
State<StateType>::findSameFockState(const StateType & ffs) {
  for (auto it = Terms.begin(); it != Terms.end(); ++it) {
    if (it->second == ffs) {
      return it;
    }
  }
  return Terms.end();
}

//+=
template<typename StateType>
State<StateType> & State<StateType>::operator+=(const StateType & rhs) {
  std::pair<std::complex<double>, StateType> toAdd(std::complex<double>(1.0, 0), rhs);
  *this += toAdd;
  return *this;
}

//+=
template<typename StateType>
State<StateType> & State<StateType>::operator+=(
    const std::pair<std::complex<double>, StateType> & rhs) {
  if (std::abs(rhs.first) < ERROR) {
    return *this;
  }
  auto it = findSameFockState(rhs.second);
  if (it == Terms.end()) {
    Terms.push_back(rhs);
  }
  else {
    it->first += rhs.first;
  }
  return *this;
}

//+=
template<typename StateType>
State<StateType> & State<StateType>::operator+=(const State<StateType> & rhs) {
  for (auto termIt = rhs.getBegin(); termIt != rhs.getEnd(); ++termIt) {
    *this += *termIt;
  }
  return *this;
}

//-=
template<typename StateType>
State<StateType> & State<StateType>::operator-=(const StateType & rhs) {
  std::pair<std::complex<double>, StateType> toAdd(std::complex<double>(-1.0, 0), rhs);
  *this += toAdd;
  return *this;
}

//-=
template<typename StateType>
State<StateType> & State<StateType>::operator-=(
    const std::pair<std::complex<double>, StateType> & rhs) {
  if (std::abs(rhs.first) < ERROR) {
    return *this;
  }
  auto it = findSameFockState(rhs.second);
  if (it == Terms.end()) {
    std::pair<std::complex<double>, StateType> copy(
        std::complex<double>(-1.0, 0) * rhs.first, rhs.second);
    Terms.push_back(copy);
  }
  else {
    it->first -= rhs.first;
  }
  return *this;
}

//-=
template<typename StateType>
State<StateType> & State<StateType>::operator-=(State<StateType> const & rhs) {
  for (auto termIt = rhs.getBegin(); termIt != rhs.getEnd(); ++termIt) {
    *this -= *termIt;
  }
  return *this;
}

//*=
template<typename StateType>
State<StateType> & State<StateType>::operator*=(std::complex<double> pref) {
  for (auto termIt = Terms.begin(); termIt != Terms.end(); ++termIt) {
    termIt->first *= pref;
  }
  return *this;
}

template<typename StateType>
bool isZeroState(std::pair<std::complex<double>, StateType> term) {
  return std::abs(term.first) < ERROR;
}

template<typename StateType>
void State<StateType>::eraseZeros() {
  Terms.erase(std::remove_if(Terms.begin(), Terms.end(), isZeroState), Terms.end());
}

//-------------------------------------------------------------------Inner Product-------

template<typename StateType>
double innerProduct(const StateType & lhs, const StateType & rhs) {
  if (lhs == rhs) {
    return 1;
  }
  return 0;
}

template<typename StateType>
std::complex<double> innerProduct(const StateType & lhs, const State<StateType> & rhs) {
  std::complex<double> ans(0, 0);
  for (auto termIt = rhs.getBegin(); termIt != rhs.getEnd(); ++termIt) {
    ans += termIt->first * innerProduct(lhs, termIt->second);
  }
  return ans;
}

template<typename StateType>
std::complex<double> innerProduct(const State<StateType> & lhs, const StateType & rhs) {
  return std::conj(innerProduct(rhs, lhs));
}

template<typename StateType>
std::complex<double> innerProduct(const State<StateType> & lhs,
                                  const State<StateType> & rhs) {
  std::complex<double> ans(0, 0);
  for (auto it = lhs.getBegin(); it != lhs.getEnd(); ++it) {
    ans += it->first * innerProduct(it->second, rhs);
  }
  return ans;
}

#endif  //QM_STATES_TEM_HPP
