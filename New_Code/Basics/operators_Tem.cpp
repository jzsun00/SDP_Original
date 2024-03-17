/*
  Jiazheng Sun
  Updated: Mar 17, 2024

  Implementations of methods in class:
  Monomial, Polynomial.
 */

#ifndef ORI_SDP_GS_OPERATORS_TEM_CPP
#define ORI_SDP_GS_OPERATORS_TEM_CPP

#include "operators.hpp"

//---------------------------------------------------------------Monomial---------------

template<typename OpType>
std::string Monomial<OpType>::toString() const {
  std::string ans;
  for (typename vector<OpType>::const_iterator it = Expr.begin(); it != Expr.end();
       ++it) {
    ans += it->toString();
  }
  return ans;
}

template<typename OpType>
Monomial<OpType> & Monomial<OpType>::operator=(Monomial<OpType> const & rhs) {
  Expr = rhs.Expr;
  return *this;
}

template<typename OpType>
Monomial<OpType> & Monomial<OpType>::operator*=(Monomial<OpType> const & rhs) {
  Expr.insert(Expr.end(), rhs.Expr.begin(), rhs.Expr.end());
  return *this;
}

template<typename OpType>
Monomial<OpType> & Monomial<OpType>::operator*=(OpType const & toAdd) {
  Expr.push_back(toAdd);
  return *this;
}

template<typename OpType>
void Monomial<OpType>::herm() {
  std::reverse(Expr.begin(), Expr.end());
  for (typename vector<OpType>::iterator it = Expr.begin(); it != Expr.end(); ++it) {
    it->herm();
  }
}

//---------------------------------------------------------------Polynomial-------------

template<typename MonomialType>
std::string Polynomial<MonomialType>::toString() const {
  std::string ans;
  for (typename vector<pair<complex<double>, MonomialType> >::const_iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    ans += "  (";
    ans += complex_toString(it->first);
    ans += ")";
    ans += it->second.toString();
    ans += "  +\n";
  }
  ans.pop_back();
  ans.pop_back();
  return ans;
}

template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator=(
    Polynomial<MonomialType> const & rhs) {
  Terms = rhs.Terms;
  return *this;
}

template<typename MonomialType>
typename vector<pair<complex<double>, MonomialType> >::iterator
Polynomial<MonomialType>::findSameMonomial(MonomialType const & mn) {
  for (typename vector<pair<complex<double>, MonomialType> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    if (it->second == mn) {
      return it;
    }
  }
  return Terms.end();
}

//+=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator+=(
    MonomialType const & rhs) {
  pair<complex<double>, MonomialType> toAdd(complex<double>(1, 0), rhs);
  *this += toAdd;
  return *this;
}

template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator+=(
    pair<complex<double>, MonomialType> const & toAdd) {
  if (std::abs(toAdd.first) < ERROR) {
    return *this;
  }
  typename vector<pair<complex<double>, MonomialType> >::iterator it =
      findSameMonomial(toAdd.second);
  if (it == Terms.end()) {
    Terms.push_back(toAdd);
  }
  else {
    it->first += toAdd.first;
  }
  return *this;
}

template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator+=(
    Polynomial<MonomialType> const & rhs) {
  for (typename vector<pair<complex<double>, MonomialType> >::const_iterator termIt =
           rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this += *termIt;
  }
  return *this;
}

//-=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator-=(
    MonomialType const & rhs) {
  pair<complex<double>, MonomialType> toAdd(complex<double>(-1, 0), rhs);
  *this += toAdd;
  return *this;
}

template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator-=(
    pair<complex<double>, MonomialType> const & toAdd) {
  if (std::abs(toAdd.first) < ERROR) {
    return *this;
  }
  pair<complex<double>, MonomialType> copy(complex<double>(-1, 0) * toAdd.first,
                                           toAdd.second);
  *this += copy;
  return *this;
}

template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator-=(
    Polynomial<MonomialType> const & rhs) {
  for (typename vector<pair<complex<double>, MonomialType> >::const_iterator termIt =
           rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this -= *termIt;
  }
  return *this;
}

//*=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator*=(
    MonomialType const & rhs) {
  pair<complex<double>, MonomialType> toAdd(complex<double>(1, 0), rhs);
  *this *= toAdd;
  return *this;
}

template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator*=(
    pair<complex<double>, MonomialType> const & rhs) {
  for (typename vector<pair<complex<double>, MonomialType> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    it->first *= rhs.first;
    it->second *= rhs.second;
  }
  return *this;
}

template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator*=(
    Polynomial<MonomialType> const & rhs) {
  Polynomial current(*this);
  for (typename vector<pair<complex<double>, MonomialType> >::const_iterator termIt =
           rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    Polynomial copy(current);
    *this += (copy *= *termIt);
  }
  *this -= current;
  return *this;
}

template<typename MonomialType>
void Polynomial<MonomialType>::herm() {
  for (typename vector<pair<complex<double>, MonomialType> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    it->first = std::conj(it->first);
    it->second.herm();
  }
}

template<typename MonomialType>
bool isZero(pair<complex<double>, MonomialType> term) {
  return std::abs(term.first) < ERROR;
}

template<typename MonomialType>
void Polynomial<MonomialType>::eraseZeros() {
  Terms.erase(std::remove_if(Terms.begin(), Terms.end(), isZero<MonomialType>),
              Terms.end());
}

#endif  //ORI_SDP_GS_OPERATORS_TEM_CPP
