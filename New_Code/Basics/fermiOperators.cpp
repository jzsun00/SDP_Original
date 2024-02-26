/*
  Implementations of methods in class:
  LadderOp, Monomial, Polynomial.
 */

#ifndef ORI_SDP_GS_FERMIOPERATORS_CPP
#define ORI_SDP_GS_FERMIOPERATORS_CPP

#include "fermiOperators.hpp"

FermiState FermiLadderOp::operator*(FermiFockState const & rhs) const {
  if (creatorF) {
    if (rhs[index] == 1) {
      return FermiState(complex<double>(0, 0), rhs);
    }
    vector<bool> ans = rhs.getNums();
    ans[index] = 1;
    size_t sum = 0;
    for (int i = 0; i < index; i++) {
      sum += rhs[i];
    }
    if (sum % 2 == 0) {
      return FermiState(complex<double>(1, 0), FermiFockState(ans));
    }
    else {
      return FermiState(complex<double>(-1, 0), FermiFockState(ans));
    }
  }
  else {
    if (rhs[index] == 0) {
      return FermiState(complex<double>(0, 0), rhs);
    }
    vector<bool> ans = rhs.getNums();
    ans[index] = 0;
    size_t sum = 0;
    for (int i = 0; i < index; i++) {
      sum += rhs[i];
    }
    if (sum % 2 == 0) {
      return FermiState(complex<double>(1, 0), FermiFockState(ans));
    }
    else {
      return FermiState(complex<double>(1, 0), FermiFockState(ans));
    }
  }
}

#endif
