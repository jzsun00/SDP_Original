/*
  Jiazheng Sun
  Updated: Jul 23, 2024

  Class Implementations:
  LadderOp<IndexType>
  SpinOp<IndexType>
  Monomial<OpType>
  Polynomial<MonomialType>
*/

#ifndef QM_OPERATORS_TEM_HPP
#define QM_OPERATORS_TEM_HPP

#include "./operators.hpp"

//----------------------------------------------------------LadderOp<IndexType>---------

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

template<typename IndexType>
LadderOp<IndexType> & LadderOp<IndexType>::operator=(const Operator<IndexType> & rhs) {
  this->index = rhs.index;
  this->isUnit = rhs.isUnit;
  this->creatorF = rhs.creatorF;
  return *this;
}

template<typename IndexType>
bool LadderOp<IndexType>::operator==(const LadderOp<IndexType> & rhs) const {
  if (this->isUnit && rhs.isUnit) {
    return true;
  }
  return (this->index == rhs.index) && (this->creatorF == rhs.creatorF);
}

//-----------------------------------------------------------SpinOp<IndexType>----------

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
SpinOp<IndexType> & SpinOp<IndexType>::operator=(const SpinOp<IndexType> & rhs) {
  this->index = rhs.index;
  this->isZ = rhs.isZ;
  this->isPlus = rhs.isPlus;
  this->isUnit = rhs.isUnit;
  return *this;
}

template<typename IndexType>
bool SpinOp<IndexType>::operator==(SpinOp<IndexType> const & rhs) const {
  return (this->index == rhs.index) && (this->isZ == rhs.isZ) &&
         (this->isPlus == rhs.isPlus) && (this->isUnit == rhs.isUnit);
}

template<typename IndexType>
void SpinOp<IndexType>::herm() {
  if (isZ || isUnit) {
    return;
  }
  else {
    isPlus ^= 1;
  }
}

//------------------------------------------------------------Monomial<OpType>----------

template<typename OpType>
std::string Monomial<OpType>::toString() const {
  std::string ans;
  for (auto it = Expr.begin(); it != Expr.end(); ++it) {
    ans += it->toString();
  }
  return ans;
}

template<typename OpType>
Monomial<OpType> & Monomial<OpType>::operator=(const Monomial<OpType> & rhs) {
  Expr = rhs.Expr;
  return *this;
}

template<typename OpType>
Monomial<OpType> & Monomial<OpType>::operator*=(const Monomial<OpType> & rhs) {
  if (rhs.getSize() == 1 && rhs[0].getIsUnit()) {
    return *this;
  }
  if (this->getSize() == 1 && this->Expr[0].getIsUnit()) {
    *this = rhs;
    return *this;
  }
  Expr.insert(Expr.end(), rhs.Expr.begin(), rhs.Expr.end());
  return *this;
}

template<typename OpType>
Monomial<OpType> & Monomial<OpType>::operator*=(const OpType & toAdd) {
  if (toAdd.getIsUnit() && Expr.size() > 0) {
    return *this;
  }
  Expr.push_back(toAdd);
  return *this;
}

template<typename OpType>
void Monomial<OpType>::herm() {
  std::reverse(Expr.begin(), Expr.end());
  for (auto it = Expr.begin(); it != Expr.end(); ++it) {
    it->herm();
  }
}

//-------------------------------------------------------Polynomial<MonomialType>-------

template<typename MonomialType>
std::string Polynomial<MonomialType>::toString() const {
  std::string ans;
  if (Terms.size() == 0) {
    return ans;
  }
  for (auto it = Terms.begin(); it != Terms.end(); ++it) {
    ans += "\t(";
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
auto Polynomial<MonomialType>::findSameMonomial(const MonomialType & mn) const {
  for (auto it = Terms.begin(); it != Terms.end(); ++it) {
    if (it->second == mn) {
      return it;
    }
  }
  return Terms.end();
}

//+=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator+=(
    const MonomialType & rhs) {
  std::pair<std::complex<double>, MonomialType> toAdd(std::complex<double>(1.0, 0), rhs);
  *this += toAdd;
  return *this;
}

//+=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator+=(
    const std::pair<std::complex<double>, MonomialType> & toAdd) {
  if (std::abs(toAdd.first) < ERROR) {
    return *this;
  }
  auto it = findSameMonomial(toAdd.second);
  if (it == Terms.end()) {
    Terms.push_back(toAdd);
  }
  else {
    it->first += toAdd.first;
  }
  return *this;
}

//+=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator+=(
    const Polynomial<MonomialType> & rhs) {
  for (auto termIt = rhs.getBegin(); termIt != rhs.getEnd(); ++termIt) {
    *this += *termIt;
  }
  return *this;
}

//-=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator-=(
    const MonomialType & rhs) {
  std::pair<std::complex<double>, MonomialType> toAdd(std::complex<double>(-1.0, 0), rhs);
  *this += toAdd;
  return *this;
}

//-=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator-=(
    const std::pair<std::complex<double>, MonomialType> & toAdd) {
  if (std::abs(toAdd.first) < ERROR) {
    return *this;
  }
  std::pair<std::complex<double>, MonomialType> copy(
      std::complex<double>(-1.0, 0) * toAdd.first, toAdd.second);
  *this += copy;
  return *this;
}

//-=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator-=(
    const Polynomial<MonomialType> & rhs) {
  for (auto termIt = rhs.getBegin(); termIt != rhs.getEnd(); ++termIt) {
    *this -= *termIt;
  }
  return *this;
}

//*=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator*=(
    const MonomialType & rhs) {
  std::pair<std::complex<double>, MonomialType> toAdd(std::complex<double>(1.0, 0), rhs);
  *this *= toAdd;
  return *this;
}

//*=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator*=(
    const std::pair<std::complex<double>, MonomialType> & rhs) {
  for (auto it = Terms.begin(); it != Terms.end(); ++it) {
    it->first *= rhs.first;
    it->second *= rhs.second;
  }
  return *this;
}

//*=
template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator*=(
    const Polynomial<MonomialType> & rhs) {
  const Polynomial<MonomialType> current(*this);
  for (auto termIt = rhs.getBegin(); termIt != rhs.getEnd(); ++termIt) {
    Polynomial<MonomialType> copy(current);
    *this += (copy *= *termIt);
  }
  *this -= current;
  return *this;
}

template<typename MonomialType>
Polynomial<MonomialType> & Polynomial<MonomialType>::operator*=(
    const std::complex<double> rhs) {
  for (auto it = Terms.begin(); it != Terms.end(); ++it) {
    it->first *= rhs;
  }
  return *this;
}

template<typename MonomialType>
void Polynomial<MonomialType>::herm() {
  for (auto it = Terms.begin(); it != Terms.end(); ++it) {
    it->first = std::conj(it->first);
    it->second.herm();
  }
}

template<typename MonomialType>
bool isZero(std::pair<std::complex<double>, MonomialType> term) {
  return std::abs(term.first) < ERROR;
}

template<typename MonomialType>
void Polynomial<MonomialType>::eraseZeros() {
  Terms.erase(std::remove_if(Terms.begin(), Terms.end(), isZero<MonomialType>),
              Terms.end());
}

#endif  //QM_OPERATORS_TEM_HPP
