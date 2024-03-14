/*
  Jiazheng Sun
  Updated: Mar 12, 2024

  Implementations of methods in class:
  FermiFockstate, FermiState, FermiBasis.
*/

#ifndef ORI_SDP_GS_SPINOPERATORS_CPP
#define ORI_SDP_GS_SPINOPERATORS_CPP

#include "spinOperators1D.hpp"

//------------------------------------------------------------------SpinZHalfOp----------

std::string SpinZHalfOp::toString() const {
  std::string ans = "SZ_{";
  ans += std::to_string(this->index);
  ans += "}";
  return ans;
}

SpinHalfState SpinZHalfOp::operator*(SpinHalfBaseState const & rhs) const {
  if (rhs[index]) {
    return SpinHalfState(complex<double>(0.5, 0), rhs);
  }
  else {
    return SpinHalfState(complex<double>(-0.5, 0), rhs);
  }
}

SpinHalfState SpinZHalfOp::operator*(SpinHalfState const & rhs) const {
  SpinHalfState ans;
  for (vector<pair<complex<double>, SpinHalfBaseState> >::const_iterator it =
           rhs.getBegin();
       it != rhs.getEnd();
       ++it) {
    ans += (((*this) * (it->second)) *= it->first);
  }
  return ans;
}

//------------------------------------------------------------------SpinUDHalfOp---------

std::string SpinUDHalfOp::toString() const {
  std::string ans = "S";
  if (this->plusF) {
    ans += "+{";
  }
  else {
    ans += "-{";
  }
  ans += std::to_string(this->index);
  ans += "}";
  return ans;
}

SpinHalfState SpinUDHalfOp::operator*(SpinHalfBaseState const & rhs) const {
  if (plusF) {
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

SpinHalfState SpinUDHalfOp::operator*(SpinHalfState const & rhs) const {
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
  for (vector<SpinHalfOp *>::const_iterator it = ExprPtr.begin(); it != ExprPtr.end();
       ++it) {
    ans = (**it) * ans;
  }
  return ans;
}

SpinHalfState SpinHalfMonomial::operator*(SpinHalfState const & rhs) const {
  SpinHalfState ans(rhs);
  for (vector<SpinHalfOp *>::const_iterator it = ExprPtr.begin(); it != ExprPtr.end();
       ++it) {
    ans = (**it) * ans;
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

vector<pair<complex<double>, SpinHalfMonomial> >::iterator
SpinHalfPolynomial::findSameMonomial(SpinHalfMonomial const & mn) {
  //std::cout << "In the dynamic dispatch area" << std::endl;
  //std::cout << "size(Terms) = " << Terms.size() << std::endl;
  //std::cout << "Terms = " << this->toString() << std::endl;
  //std::cout << "Trying to find same monomial for " << mn.toString() << std::endl;
  /*
  for (typename vector<pair<complex<double>, SpinHalfMonomial> >::iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    if (it->second == mn) {
      return it;
    }
  }
  */
  //std::cout << "Returning Terms.end()" << std::endl;
  return Terms.end();
}

/*
SpinHalfPolynomial & SpinHalfPolynomial::operator+=(
    pair<complex<double>, SpinHalfMonomial> const & toAdd) {
  std::cout << "In the dynamic dispatch +=" << std::endl;
  if (std::abs(toAdd.first) < ERROR) {
    return *this;
  }
  typename vector<pair<complex<double>, SpinHalfMonomial> >::iterator it =
      findSameMonomial(toAdd.second);
  if (it == Terms.end()) {
    Terms.push_back(toAdd);
  }
  else {
    it->first += toAdd.first;
  }
  return *this;
}
*/
#endif
