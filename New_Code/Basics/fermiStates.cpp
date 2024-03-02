/*
  Implementations of methods in class:
  LadderOp, Monomial, Polynomial.
 */

#ifndef ORI_SDP_GS_FERMISTATES_CPP
#define ORI_SDP_GS_FERMISTATES_CPP

#include "fermiStates.hpp"

//-------------------------------------------------------------FermiFockState-----------

using std::vector;

std::string FermiFockState::toString() {
  std::string ans = "|";
  for (vector<bool>::iterator it = Nums.begin(); it != Nums.end(); ++it) {
    ans += " ";
    ans += std::to_string(*it);
    ans += ",";
  }
  ans.pop_back();
  ans += " >";
  return ans;
}

FermiFockState & FermiFockState::operator=(FermiFockState const & rhs) {
  Nums = rhs.Nums;
  return *this;
}

//--------------------------------------------------------------FermiState--------------

std::string FermiState::toString() {
  std::string ans = "";
  for (vector<pair<complex<double>, FermiFockState> >::iterator it = Terms.begin();
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

pair<complex<double>, FermiFockState> FermiState::operator[](size_t n) const {
  return Terms.at(n);
}

FermiState & FermiState::operator=(FermiState const & rhs) {
  Terms = rhs.Terms;
  return *this;
}

vector<pair<complex<double>, FermiFockState> >::iterator FermiState::findSameFockState(
    FermiFockState const & ffs) {
  for (vector<pair<complex<double>, FermiFockState> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    if (it->second == ffs) {
      return it;
    }
  }
  return Terms.end();
}

FermiState & FermiState::operator+=(FermiFockState const & rhs) {
  pair<complex<double>, FermiFockState> toAdd(complex<double>(1, 0), rhs);
  *this += toAdd;
  return *this;
}

FermiState & FermiState::operator+=(pair<complex<double>, FermiFockState> const & rhs) {
  if (std::abs(rhs.first) < std::pow(10, -12)) {
    return *this;
  }
  vector<pair<complex<double>, FermiFockState> >::iterator it =
      findSameFockState(rhs.second);
  if (it == Terms.end()) {
    Terms.push_back(rhs);
  }
  else {
    it->first += rhs.first;
  }
  return *this;
}

FermiState & FermiState::operator+=(FermiState & rhs) {
  for (vector<pair<complex<double>, FermiFockState> >::iterator termIt = rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this += *termIt;
  }
  return *this;
}

FermiState & FermiState::operator-=(FermiFockState const & rhs) {
  pair<complex<double>, FermiFockState> toAdd(complex<double>(1, 0), rhs);
  *this -= toAdd;
  return *this;
}

FermiState & FermiState::operator-=(pair<complex<double>, FermiFockState> const & rhs) {
  if (std::abs(rhs.first) < std::pow(10, -12)) {
    return *this;
  }
  vector<pair<complex<double>, FermiFockState> >::iterator it =
      findSameFockState(rhs.second);
  if (it == Terms.end()) {
    pair<complex<double>, FermiFockState> copy(complex<double>(-1, 0) * rhs.first,
                                               rhs.second);
    Terms.push_back(copy);
  }
  else {
    it->first -= rhs.first;
  }
  return *this;
}

FermiState & FermiState::operator-=(FermiState & rhs) {
  for (vector<pair<complex<double>, FermiFockState> >::iterator termIt = rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this -= *termIt;
  }
  return *this;
}

FermiState & FermiState::operator*=(complex<double> pref) {
  for (vector<pair<complex<double>, FermiFockState> >::iterator termIt = Terms.begin();
       termIt != Terms.end();
       ++termIt) {
    termIt->first *= pref;
  }
  return *this;
}

bool isZero(pair<complex<double>, FermiFockState> term) {
  return std::abs(term.first) < std::pow(10, -12);
}

void FermiState::eraseZeros() {
  Terms.erase(std::remove_if(Terms.begin(), Terms.end(), isZero), Terms.end());
}

//-------------------------------------------------------------FermiBasis---------------

void FermiBasis::init() {
  size_t total = std::pow(2, Sites);
  for (size_t i = 0; i < total; i++) {
    std::bitset<32> bits(i);
    //std::cout << "i = " << i << ", bitset = " << bits << std::endl;
    vector<bool> newState(Sites);
    for (size_t j = 0; j < Sites; j++) {
      newState[Sites - j - 1] = bits[j];
      //std::cout << "j + 32 - Sites = " << j + 32 - Sites << std::endl;
      //std::cout << "bits[j + 32 - Sites] = " << bits[j + 32 - Sites] << std::endl;
    }
    FermiFockState toAdd(newState);
    //std::cout << "toAdd = " << toAdd.toString() << std::endl;
    States.push_back(toAdd);
  }
}

std::string FermiBasis::toString() {
  std::string ans = "Sites = ";
  ans += std::to_string(Sites);
  ans += "\nFull Basis:\n";
  for (vector<FermiFockState>::iterator it = States.begin(); it != States.end(); ++it) {
    ans += it->toString();
    ans += "\n";
  }
  return ans;
}

//-------------------------------------------------------------Other Functions---------------

double innerProduct(FermiFockState lhs, FermiFockState rhs) {
  if (lhs == rhs) {
    return 1;
  }
  return 0;
}

complex<double> innerProduct(FermiFockState lhs, FermiState rhs) {
  complex<double> ans(0, 0);
  for (vector<pair<complex<double>, FermiFockState> >::iterator termIt = rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    ans += termIt->first * innerProduct(lhs, termIt->second);
  }
  return ans;
}

complex<double> innerProduct(FermiState lhs, FermiFockState rhs) {
  return std::conj(innerProduct(rhs, lhs));
}

#endif
