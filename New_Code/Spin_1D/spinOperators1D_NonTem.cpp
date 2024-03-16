/*
  Jiazheng Sun
  Updated: Mar 16, 2024

  Implementations of methods in class:
  SpinHalfOp, SpinHalfMonomial, SpinHalfPolynomial.
*/

#ifndef ORI_SDP_GS_SPINOPERATORS1D_NONTEM_CPP
#define ORI_SDP_GS_SPINOPERATORS1D_NONTEM_CPP

#include "spinOperators1D.hpp"

//------------------------------------------------------------------SpinHalfOp----------

SpinHalfState SpinHalfOp::operator*(SpinHalfBaseState const & rhs) const {
  if (isZ) {
    if (rhs[index]) {
      return SpinHalfState(complex<double>(0.5, 0), rhs);
    }
    else {
      return SpinHalfState(complex<double>(-0.5, 0), rhs);
    }
  }
  else {
    if (isPlus) {
      if (rhs[index]) {
        return SpinHalfState(complex<double>(0, 0), rhs);
      }
      vector<bool> Nums = rhs.getNums();
      Nums[index] = true;
      return SpinHalfState(complex<double>(1.0, 0), SpinHalfBaseState(Nums));
    }
    else {
      if (!rhs[index]) {
        return SpinHalfState(complex<double>(0, 0), rhs);
      }
      vector<bool> Nums = rhs.getNums();
      Nums[index] = false;
      return SpinHalfState(complex<double>(1.0, 0), SpinHalfBaseState(Nums));
    }
  }
}

SpinHalfState SpinHalfOp::operator*(SpinHalfState const & rhs) const {
  SpinHalfState ans;
  for (vector<pair<complex<double>, SpinHalfBaseState> >::const_iterator it =
           rhs.getBegin();
       it != rhs.getEnd();
       ++it) {
    ans += (((*this) * (it->second)) *= it->first);
  }
  return ans;
}

//-------------------------------------------------------------------SpinMonomial--------

SpinHalfState SpinHalfMonomial::operator*(SpinHalfBaseState const & rhs) const {
  SpinHalfState ans(rhs);
  for (vector<SpinHalfOp>::const_iterator it = Expr.begin(); it != Expr.end(); ++it) {
    ans = (*it) * ans;
  }
  return ans;
}

SpinHalfState SpinHalfMonomial::operator*(SpinHalfState const & rhs) const {
  SpinHalfState ans(rhs);
  for (vector<SpinHalfOp>::const_iterator it = Expr.begin(); it != Expr.end(); ++it) {
    ans = (*it) * ans;
  }
  return ans;
}

//------------------------------------------------------------------SpinPolynomial-------

SpinHalfState SpinHalfPolynomial::operator*(SpinHalfBaseState const & rhs) const {
  SpinHalfState ans;
  for (vector<pair<complex<double>, SpinHalfMonomial> >::const_iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    ans += (((it->second) * rhs) *= it->first);
  }
  return ans;
}

SpinHalfState SpinHalfPolynomial::operator*(SpinHalfState const & rhs) const {
  SpinHalfState ans;
  for (vector<pair<complex<double>, SpinHalfMonomial> >::const_iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    ans += (((it->second) * rhs) *= it->first);
  }
  return ans;
}

#endif  //ORI_SDP_GS_SPINOPERATORS1D_NONTEM_CPP
