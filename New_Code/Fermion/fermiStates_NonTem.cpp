/*
  Jiazheng Sun
  Updated: Jul 31, 2024
  
  Class Implementations:
  FermiFockstate
  FermiState
  FermiBasis.
*/

#ifndef QM_FERMI_STATES_NONTEM_CPP
#define QM_FERMI_STATES_NONTEM_CPP

#include "./fermiStates.hpp"

using std::vector;

//-------------------------------------------------------------FermiBasis----------------

void Fermi1DBasis::init() {
  size_t total = std::pow(2, Sites);
  for (size_t i = 0; i < total; i++) {
    std::bitset<32> bits(i);
    vector<bool> newState(Sites);
    for (size_t j = 0; j < Sites; j++) {
      newState[Sites - j - 1] = bits[j];
    }
    FermiFockState toAdd(newState);
    States.push_back(toAdd);
  }
}

std::string Fermi1DBasis::toString() {
  std::string ans = "1D Fermi System Basis\nSites = ";
  ans += std::to_string(Sites);
  ans += "\nFull Basis:\n";
  for (vector<FermiFockState>::const_iterator it = States.begin(); it != States.end();
       ++it) {
    ans += it->toString();
    ans += "\n";
  }
  return ans;
}

#endif  //QM_FERMI_STATES_NONTEM_CPP
