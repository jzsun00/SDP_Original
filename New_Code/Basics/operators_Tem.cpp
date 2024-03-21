/*
  Jiazheng Sun
  Updated: Mar 18, 2024

  Implementations of methods in class:
  LadderOp, SpinOp, Monomial, Polynomial.
 */

#ifndef ORI_SDP_GS_OPERATORS_TEM_CPP
#define ORI_SDP_GS_OPERATORS_TEM_CPP

#include "./operators.hpp"

//---------------------------------------------------------------LadderOp---------------

template<typename IndexType>
std::string LadderOp<IndexType>::toString() const {
  if (isUnit) {
    std::string ans = "e";
    return ans;
  }
  std::string ans = "a_{";
  ans += indexToString();
  ans += "}";
  if (creatorF) {
    ans += "{+}";
  }
  return ans;
}

//------------------------------------------------------------------SpinOp--------------

template<typename IndexType>
std::string SpinOp<IndexType>::toString() const {
  std::string ans = "S";
  if (isZ) {
    ans += "z";
  }
  else {
    if (isPlus) {
      ans += "+";
    }
    else {
      ans += "-";
    }
  }
  ans += "_{";
  ans += indexToString();
  ans += "}";
  return ans;
}

template<typename IndexType>
SpinOp<IndexType> & SpinOp<IndexType>::operator=(SpinOp<IndexType> const & rhs) {
  this->index = rhs.index;
  this->isZ = rhs.isZ;
  this->isPlus = rhs.isPlus;
  return *this;
}

template<typename IndexType>
bool SpinOp<IndexType>::operator==(SpinOp<IndexType> const & rhs) const {
  return (this->index == rhs.index) && (this->isZ == rhs.isZ) &&
         (this->isPlus == rhs.isPlus);
}

template<typename IndexType>
void SpinOp<IndexType>::herm() {
  if (isZ) {
    return;
  }
  else {
    isPlus ^= 1;
  }
}

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
  if (rhs.getSize() == 1 && rhs[0].getIsUnit()) {
    return *this;
  }
  Expr.insert(Expr.end(), rhs.Expr.begin(), rhs.Expr.end());
  return *this;
}

template<typename OpType>
Monomial<OpType> & Monomial<OpType>::operator*=(OpType const & toAdd) {
  if (toAdd.getIsUnit() && Expr.size() > 0) {
    return *this;
  }
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
  if (Terms.size() == 0) {
    return ans;
  }
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
Polynomial<MonomialType> & Polynomial<MonomialType>::operator*=(complex<double> rhs) {
  for (typename vector<pair<complex<double>, MonomialType> >::iterator it = Terms.begin();
       it != Terms.end();
       ++it) {
    it->first *= rhs;
  }
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
