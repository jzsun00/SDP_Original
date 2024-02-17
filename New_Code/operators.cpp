/*
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

std::string Monomial::toString() {
  std::string ans = "";
  for (vector<LadderOp>::iterator it = Expr.begin(); it != Expr.end(); ++it) {
    ans += it->toString();
  }
  return ans;
}

Monomial & Monomial::operator=(Monomial const & rhs) {
  Expr = rhs.Expr;
  return *this;
}

Monomial & Monomial::operator*=(Monomial const & rhs) {
  Expr.insert(Expr.end(), rhs.Expr.begin(), rhs.Expr.end());
  return *this;
}
Monomial & Monomial::operator*=(LadderOp const & toAdd) {
  Expr.push_back(toAdd);
  return *this;
}
void Monomial::herm() {
  std::reverse(Expr.begin(), Expr.end());
  for (vector<LadderOp>::iterator it = Expr.begin(); it != Expr.end(); ++it) {
    it->herm();
  }
}

//---------------------------------------------------------------Polynomial-------------

std::string Polynomial::toString() {
  std::string ans = "";
  for (vector<pair<complex<double>, Monomial> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    ans += "  (";
    ans += std::to_string(it->first.real());
    ans += " + ";
    ans += std::to_string(it->first.imag());
    ans += ")";
    ans += it->second.toString();
    ans += "  +";
  }
  ans.pop_back();
  return ans;
}

vector<pair<complex<double>, Monomial> >::iterator Polynomial::findSameMonomial(
    Monomial const & mn) {
  for (vector<pair<complex<double>, Monomial> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    if (it->second == mn) {
      return it;
    }
  }
  return Terms.end();
}

Polynomial & Polynomial::operator+=(pair<complex<double>, Monomial> const & toAdd) {
  vector<pair<complex<double>, Monomial> >::iterator it = findSameMonomial(toAdd.second);
  if (it == Terms.end()) {
    Terms.push_back(toAdd);
  }
  else {
    it->first += toAdd.first;
  }
  return *this;
}

#endif
