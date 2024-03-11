/*
  Jiazheng Sun
  Updated: Mar 9, 2024
  Implementations of methods in class:
  LadderOp, Monomial, Polynomial.
 */

#ifndef ORI_SDP_GS_OPERATORS_CPP
#define ORI_SDP_GS_OPERATORS_CPP

#include "operators.hpp"

//---------------------------------------------------------------LadderOp---------------

std::string LadderOp::toString() const {
  std::string ans = "a_{";
  ans += std::to_string(index);
  ans += "}";
  if (creatorF) {
    ans += "{+}";
  }
  return ans;
}

LadderOp & LadderOp::operator=(LadderOp const & rhs) {
  index = rhs.index;
  creatorF = rhs.creatorF;
  return *this;
}

bool LadderOp::operator==(LadderOp const & rhs) const {
  return creatorF == rhs.creatorF && index == rhs.index;
}

bool LadderOp::operator<(LadderOp const & rhs) const {
  if (creatorF != rhs.creatorF) {
    throw std::invalid_argument("Two operators are not the same kind!\n");
  }
  else {
    return index < rhs.index;
  }
}

//---------------------------------------------------------------Monomial---------------

template<typename LadderType>
std::string Monomial<LadderType>::toString() {
  std::string ans = "";
  for (typename vector<LadderType>::iterator it = Expr.begin(); it != Expr.end(); ++it) {
    ans += it->toString();
  }
  return ans;
}

template<typename LadderType>
Monomial<LadderType> & Monomial<LadderType>::operator=(Monomial<LadderType> const & rhs) {
  Expr = rhs.Expr;
  return *this;
}

template<typename LadderType>
Monomial<LadderType> & Monomial<LadderType>::operator*=(
    Monomial<LadderType> const & rhs) {
  Expr.insert(Expr.end(), rhs.Expr.begin(), rhs.Expr.end());
  return *this;
}

template<typename LadderType>
Monomial<LadderType> & Monomial<LadderType>::operator*=(LadderType const & toAdd) {
  Expr.push_back(toAdd);
  return *this;
}

template<typename LadderType>
void Monomial<LadderType>::herm() {
  std::reverse(Expr.begin(), Expr.end());
  for (typename vector<LadderType>::iterator it = Expr.begin(); it != Expr.end(); ++it) {
    it->herm();
  }
}

//---------------------------------------------------------------Polynomial-------------

template<typename LadderType>
std::string Polynomial<LadderType>::toString() {
  std::string ans = "";
  for (typename vector<pair<complex<double>, Monomial<LadderType> > >::iterator it =
           Terms.begin();
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

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator=(
    Polynomial<LadderType> const & rhs) {
  Terms = rhs.Terms;
  return *this;
}

template<typename LadderType>
typename vector<pair<complex<double>, Monomial<LadderType> > >::iterator
Polynomial<LadderType>::findSameMonomial(Monomial<LadderType> const & mn) {
  for (typename vector<pair<complex<double>, Monomial<LadderType> > >::iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    if (it->second == mn) {
      return it;
    }
  }
  return Terms.end();
}

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator+=(
    Monomial<LadderType> const & rhs) {
  pair<complex<double>, Monomial<LadderType> > toAdd(complex<double>(1, 0), rhs);
  *this += toAdd;
  return *this;
}

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator+=(
    pair<complex<double>, Monomial<LadderType> > const & toAdd) {
  if (std::abs(toAdd.first) < std::pow(10, -12)) {
    return *this;
  }
  typename vector<pair<complex<double>, Monomial<LadderType> > >::iterator it =
      findSameMonomial(toAdd.second);
  if (it == Terms.end()) {
    Terms.push_back(toAdd);
  }
  else {
    it->first += toAdd.first;
  }
  return *this;
}

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator+=(
    Polynomial<LadderType> const & rhs) {
  for (typename vector<pair<complex<double>, Monomial<LadderType> > >::iterator termIt =
           rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this += *termIt;
  }
  return *this;
}

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator-=(
    Monomial<LadderType> const & rhs) {
  pair<complex<double>, Monomial<LadderType> > toAdd(complex<double>(-1, 0), rhs);
  *this += toAdd;
  return *this;
}

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator-=(
    pair<complex<double>, Monomial<LadderType> > const & toAdd) {
  if (std::abs(toAdd.first) < std::pow(10, -12)) {
    return *this;
  }
  pair<complex<double>, Monomial<LadderType> > copy(-1 * toAdd.first, toAdd.sencond);
  *this += copy;
  return *this;
}

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator-=(
    Polynomial<LadderType> const & rhs) {
  for (typename vector<pair<complex<double>, Monomial<LadderType> > >::iterator termIt =
           rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this -= *termIt;
  }
  return *this;
}

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator*=(
    Monomial<LadderType> const & rhs) {
  pair<complex<double>, Monomial<LadderType> > toAdd(complex<double>(1, 0), rhs);
  *this *= toAdd;
  return *this;
}

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator*=(
    pair<complex<double>, Monomial<LadderType> > const & rhs) {
  for (typename vector<pair<complex<double>, Monomial<LadderType> > >::iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    it->first *= rhs.first;
    it->second *= rhs.second;
  }
  return *this;
}

template<typename LadderType>
Polynomial<LadderType> & Polynomial<LadderType>::operator*=(
    Polynomial<LadderType> const & rhs) {
  Polynomial current(*this);
  for (typename vector<pair<complex<double>, Monomial<LadderType> > >::iterator termIt =
           rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    Polynomial copy(current);
    *this += (copy *= *termIt);
  }
  *this -= current;
  return *this;
}

template<typename LadderType>
void Polynomial<LadderType>::herm() {
  for (typename vector<pair<complex<double>, Monomial<LadderType> > >::iterator it =
           Terms.begin();
       it != Terms.end();
       ++it) {
    it->first = std::conj(it->first);
    it->second.herm();
  }
}

template<typename LadderType>
bool isZero(pair<complex<double>, Monomial<LadderType> > term) {
  return std::abs(term.first) < std::pow(10, -12);
}

template<typename LadderType>
void Polynomial<LadderType>::eraseZeros() {
  /*
  for (vector<pair<complex<double>, Monomial> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    std::cout << "std::abs(it->first) = " << std::abs(it->first) << std::endl;
    if (std::abs(it->first) < std::pow(10, -9)) {
      Terms.erase(it);
    }
  }
  */
  Terms.erase(std::remove_if(Terms.begin(), Terms.end(), isZero), Terms.end());
}

#endif  //ORI_SDP_GS_OPERATORS_CPP
