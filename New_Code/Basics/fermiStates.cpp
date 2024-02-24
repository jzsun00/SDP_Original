/*
  Implementations of methods in class:
  LadderOp, Monomial, Polynomial.
 */

#ifndef ORI_SDP_GS_FERMISTATES_CPP
#define ORI_SDP_GS_FERMISTATES_CPP

#include "fermiStates.hpp"

double innerProduct(FermiFockState lhs, FermiFockState rhs) {
  if (lhs == rhs) {
    return 1;
  }
  return 0;
}

complex<double> innerProduct(FermiFockState lhs, FermiState rhs) {
  complex<double> ans(0, 0);
  for (vector<pair<complex<double>, FermiFockState> >::iterator termIt = rhs.getBegin();
  termIt != rhs.getEnd(); ++termIt) {
    ans += termIt->first * innerProduct(lhs, termIt->second);
  }
  return ans;
}

complex<double> innerProduct(FermiState lhs, FermiFockState rhs) {
  return std::conj(innerProduct(rhs, lhs));
}




#endif
