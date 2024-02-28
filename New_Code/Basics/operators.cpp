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
    ans += "  +\n";
  }
  ans.pop_back();
  ans.pop_back();
  return ans;
}

Polynomial & Polynomial::operator=(Polynomial const & rhs) {
  Terms = rhs.Terms;
  return *this;
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

Polynomial & Polynomial::operator+=(Monomial const & rhs) {
  pair<complex<double>, Monomial> toAdd(complex<double>(1, 0), rhs);
  *this += toAdd;
  return *this;
}

Polynomial & Polynomial::operator+=(pair<complex<double>, Monomial> const & toAdd) {
  if (std::abs(toAdd.first) < std::pow(10, -12)) {
    return *this;
  }
  vector<pair<complex<double>, Monomial> >::iterator it = findSameMonomial(toAdd.second);
  if (it == Terms.end()) {
    Terms.push_back(toAdd);
  }
  else {
    it->first += toAdd.first;
  }
  return *this;
}

Polynomial & Polynomial::operator+=(Polynomial & rhs) {
  for (vector<pair<complex<double>, Monomial> >::iterator termIt = rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this += *termIt;
  }
  return *this;
}

Polynomial & Polynomial::operator-=(Monomial const & rhs) {
  pair<complex<double>, Monomial> toAdd(complex<double>(1, 0), rhs);
  *this -= toAdd;
  return *this;
}

Polynomial & Polynomial::operator-=(pair<complex<double>, Monomial> const & toAdd) {
  if (std::abs(toAdd.first) < std::pow(10, -12)) {
    return *this;
  }
  vector<pair<complex<double>, Monomial> >::iterator it = findSameMonomial(toAdd.second);
  if (it == Terms.end()) {
    pair<complex<double>, Monomial> copy(complex<double>(-1, 0) * toAdd.first,
                                         toAdd.second);
    Terms.push_back(copy);
  }
  else {
    it->first -= toAdd.first;
  }
  return *this;
}

Polynomial & Polynomial::operator-=(Polynomial & rhs) {
  for (vector<pair<complex<double>, Monomial> >::iterator termIt = rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    *this -= *termIt;
  }
  return *this;
}

Polynomial & Polynomial::operator*=(Monomial const & rhs) {
  pair<complex<double>, Monomial> toAdd(complex<double>(1, 0), rhs);
  for (vector<pair<complex<double>, Monomial> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    it->first *= toAdd.first;
    it->second *= toAdd.second;
  }
  return *this;
}

Polynomial & Polynomial::operator*=(pair<complex<double>, Monomial> const & rhs) {
  for (vector<pair<complex<double>, Monomial> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    it->first *= rhs.first;
    it->second *= rhs.second;
  }
  return *this;
}

Polynomial & Polynomial::operator*=(Polynomial & rhs) {
  Polynomial current(*this);
  for (vector<pair<complex<double>, Monomial> >::iterator termIt = rhs.getBegin();
       termIt != rhs.getEnd();
       ++termIt) {
    Polynomial copy(current);
    *this += (copy *= *termIt);
  }
  *this -= current;
  return *this;
}

void Polynomial::herm() {
  for (vector<pair<complex<double>, Monomial> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    it->first = std::conj(it->first);
    it->second.herm();
  }
}

bool isZero(pair<complex<double>, Monomial> term) {
  return std::abs(term.first) < std::pow(10, -12);
}

void Polynomial::eraseZeros() {
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

#endif
