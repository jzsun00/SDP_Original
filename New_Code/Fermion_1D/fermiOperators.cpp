/*
  Jiazheng Sun
  Updated: Mar 9, 2024

  Implementations of methods in class:
  FermiLadderOp, FermiMonomial, FermiPolynomial.
 */

#ifndef ORI_SDP_GS_FERMIOPERATORS_CPP
#define ORI_SDP_GS_FERMIOPERATORS_CPP

#include "fermiOperators.hpp"

//------------------------------------------------------------------FermiLadderOp--------

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

//------------------------------------------------------------------FermiMonomial--------
/*
FermiState FermiMonomial::operator*(FermiFockState const & rhs) const {
  FermiState ans(rhs);
  for (int i = Expr.size() - 1; i >= 0; i--) {
    FermiLadderOp op(Expr[i]);
    ans = op * ans;
  }
  return ans;
}
*/

//-----------------------------------------------------------------FermiPolynomial-------
/*
FermiState FermiPolynomial::operator*(FermiFockState const & rhs) const {
  FermiState ans;
  for (size_t i = 0; i < Terms.size(); i++) {
    FermiMonomial mn(Terms[i].second);
    FermiState newTerm = mn * rhs;
    newTerm *= Terms[i].first;
    ans += newTerm;
  }
  return ans;
}
*/

#endif  //ORI_SDP_GS_FERMIOPERATORS_CPP
