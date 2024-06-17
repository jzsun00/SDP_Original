/*
  Jiazheng Sun
  Updated: Jun 17, 2024
  
  Class Implementations:
  SpinHalfBaseState1D
  SpinHalfState1D
  SpinHalfBasis1D
*/

#ifndef QM_SPINSTATES1D_NONTEM_CPP
#define QM_SPINSTATES1D_NONTEM_CPP

#include <cstddef>

#include "./spinStates1D.hpp"

//-------------------------------------------------------------SpinHalfState1D-----------

SpinHalfState1D & SpinHalfState1D::operator=(const SpinHalfState1D & rhs) {
  Terms = rhs.Terms;
  return *this;
}

//-------------------------------------------------------------SpinHalfBasis1D-----------

void SpinHalfBasis1D::init() {
  size_t total = std::pow(2, SitesNum);
  for (size_t i = 0; i < total; i++) {
    std::bitset<32> bits(i);
    vector<bool> newState(SitesNum);
    for (size_t j = 0; j < SitesNum; j++) {
      newState[SitesNum - j - 1] = bits[j];
    }
    SpinHalfBaseState1D toAdd(newState);
    States.push_back(toAdd);
  }
}

void SpinHalfBasis1D::init(int SzTotal) {
  size_t total = std::pow(2, SitesNum);
  for (size_t i = 0; i < total; i++) {
    std::bitset<32> bits(i);
    int Sz = 0;
    for (size_t j = 0; j < SitesNum; j++) {
      if (bits[j]) {
        Sz += 1;
      }
      else {
        Sz -= 1;
      }
    }
    if (Sz != SzTotal) {
      continue;
    }
    vector<bool> newState(SitesNum);
    for (size_t j = 0; j < SitesNum; j++) {
      newState[SitesNum - j - 1] = bits[j];
    }

    SpinHalfBaseState1D toAdd(newState);
    States.push_back(toAdd);
  }
}

std::string SpinHalfBasis1D::toString() {
  std::string ans = "SitesNum = ";
  ans += std::to_string(SitesNum);
  ans += "\nFull Basis:\n";
  for (vector<SpinHalfBaseState1D>::const_iterator it = States.begin();
       it != States.end();
       ++it) {
    ans += it->toString();
    ans += "\n";
  }
  return ans;
}

#endif  //QM_SPINSTATES1D_NONTEM_CPP
