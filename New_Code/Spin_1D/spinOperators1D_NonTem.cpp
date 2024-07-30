/*
  Jiazheng Sun
  Updated: Jul 30, 2024

  Class Implementations:
  SpinHalfOp1D
  SpinHalfMonomial1D
  SpinHalfPolynomial1D
*/

#ifndef QM_SPIN_OPERATORS_1D_NONTEM_CPP
#define QM_SPIN_OPERATORS_1D_NONTEM_CPP

#include "spinOperators1D.hpp"

using std::complex;
using std::pair;
using std::vector;

//--------------------------------------------------------------SpinHalfOp1D-------------

SpinHalfState1D SpinHalfOp1D::operator*(SpinHalfBaseState1D const & rhs) const {
  if (isZ) {
    if (rhs[index]) {
      return SpinHalfState1D(complex<double>(0.5, 0), rhs);
    }
    else {
      return SpinHalfState1D(complex<double>(-0.5, 0), rhs);
    }
  }
  else {
    if (isPlus) {
      if (rhs[index]) {
        return SpinHalfState1D(complex<double>(0, 0), rhs);
      }
      vector<bool> Nums = rhs.getAllNums();
      Nums[index] = true;
      return SpinHalfState1D(complex<double>(1.0, 0), SpinHalfBaseState1D(Nums));
    }
    else {
      if (!rhs[index]) {
        return SpinHalfState1D(complex<double>(0, 0), rhs);
      }
      vector<bool> Nums = rhs.getAllNums();
      Nums[index] = false;
      return SpinHalfState1D(complex<double>(1.0, 0), SpinHalfBaseState1D(Nums));
    }
  }
}

SpinHalfState1D SpinHalfOp1D::operator*(SpinHalfState1D const & rhs) const {
  SpinHalfState1D ans;
  for (vector<pair<complex<double>, SpinHalfBaseState1D> >::const_iterator it =
           rhs.getBegin();
       it != rhs.getEnd();
       ++it) {
    ans += (((*this) * (it->second)) *= it->first);
  }
  return ans;
}

//-----------------------------------------------------------SpinHalfMonomial1D----------

SpinHalfState1D SpinHalfMonomial1D::operator*(SpinHalfBaseState1D const & rhs) const {
  SpinHalfState1D ans(rhs);
  for (vector<SpinHalfOp1D>::const_iterator it = Expr.begin(); it != Expr.end(); ++it) {
    ans = (*it) * ans;
  }
  return ans;
}

SpinHalfState1D SpinHalfMonomial1D::operator*(SpinHalfState1D const & rhs) const {
  SpinHalfState1D ans(rhs);
  for (vector<SpinHalfOp1D>::const_iterator it = Expr.begin(); it != Expr.end(); ++it) {
    ans = (*it) * ans;
  }
  return ans;
}

//----------------------------------------------------------SpinHalfPolynomial1D---------

SpinHalfPolynomial1D & SpinHalfPolynomial1D::operator=(SpinHalfPolynomial1D const & rhs) {
  Terms = rhs.Terms;
  return *this;
}

SpinHalfState1D SpinHalfPolynomial1D::operator*(SpinHalfBaseState1D const & rhs) const {
  SpinHalfState1D ans;
  for (vector<pair<complex<double>, SpinHalfMonomial1D> >::const_iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    ans += (((it->second) * rhs) *= it->first);
  }
  return ans;
}

SpinHalfState1D SpinHalfPolynomial1D::operator*(SpinHalfState1D const & rhs) const {
  SpinHalfState1D ans;
  for (vector<pair<complex<double>, SpinHalfMonomial1D> >::const_iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    ans += (((it->second) * rhs) *= it->first);
  }
  return ans;
}

void SpinHalfPolynomial1D::operate(const SpinHalfState1D & rhs,
                                   SpinHalfState1D & result) const {
  result.clear();
  for (vector<pair<complex<double>, SpinHalfMonomial1D> >::const_iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    result += (((it->second) * rhs) *= it->first);
  }
}

#endif  //QM_SPIN_OPERATORS_1D_NONTEM_CPP
