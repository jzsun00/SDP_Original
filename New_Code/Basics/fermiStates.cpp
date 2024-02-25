/*
  Implementations of methods in class:
  LadderOp, Monomial, Polynomial.
 */

#ifndef ORI_SDP_GS_FERMISTATES_CPP
#define ORI_SDP_GS_FERMISTATES_CPP

#include "fermiStates.hpp"

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

bool FermiFockState::operator==(FermiFockState const & rhs) const {
  return Nums == rhs.Nums;
}

//-------------------------------------------

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

//------------------------------------------------

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
