/*
  Jiazheng Sun
  Updated: Jul 31, 2024
  
  Class Implementations:
  FermiLadderOp<IndexType>
  FermiMonomial<OpType>
  FermiPolynomial<MonomialType>
*/

#ifndef QM_FERMI_OPERATORS_TEM_HPP
#define QM_FERMI_OPERATORS_TEM_HPP

#include <cstddef>

#include "fermiOperators.hpp"

//------------------------------------------------------FermiLadderOp<IndexType>---------

template<typename IndexType>
FermiLadderOp<IndexType> & FermiLadderOp<IndexType>::operator=(
    const FermiLadderOp<IndexType> & rhs) {
  if (this != &rhs) {
    this->index = rhs.index;
    this->isUnit = rhs.isUnit;
    this->creatorF = rhs.creatorF;
  }
  return *this;
}

template<typename IndexType>
bool FermiLadderOp<IndexType>::operator>(const FermiLadderOp & rhs) const {
  return !(*this < rhs || *this == rhs);
}

//---------------------------------------------------------FermiMonomial<OpType>---------

template<typename OpType>
FermiMonomial<OpType> & FermiMonomial<OpType>::operator=(
    const FermiMonomial<OpType> & rhs) {
  if (this != &rhs) {
    this->Expr = rhs.Expr;
  }
  return *this;
}

template<typename OpType>
void FermiMonomial<OpType>::reverse() {
  std::reverse(this->Expr.begin(), this->Expr.end());
}

template<typename OpType>
int FermiMonomial<OpType>::findWrongOrder() const {
  for (size_t i = 0; i < this->Expr.size() - 1; i++) {
    if (this->Expr[i] > this->Expr[i + 1]) {
      return i;
    }
  }
  return -1;
}

template<typename OpType>
bool FermiMonomial<OpType>::isNorm() const {
  if (this->Expr.size() <= 1) {
    return true;
  }
  if (findWrongOrder() < 0) {
    return true;
  }
  return false;
}

template<typename OpType>
FermiMonomial<OpType> FermiMonomial<OpType>::sliceExprStart(size_t index) const {
  return FermiMonomial<OpType>(
      std::vector<OpType>(this->Expr.begin(), this->Expr.begin() + index));
}

template<typename OpType>
FermiMonomial<OpType> FermiMonomial<OpType>::sliceExprEnd(size_t index) const {
  return FermiMonomial<OpType>(
      std::vector<OpType>(this->Expr.begin() + index, this->Expr.end()));
}

template<typename OpType>
void FermiMonomial<OpType>::moveIndex(int offset) {
  size_t len = this->getSize();
  for (size_t i = 0; i < len; i++) {
    this->Expr[i].moveIndex(offset);
  }
}

template<typename OpType>
bool FermiMonomial<OpType>::equiv(const FermiMonomial<OpType> & rhs) const {
  if (this->getSize() != rhs.getSize()) {
    return false;
  }
  FermiMonomial<OpType> rhsCopy(rhs);
  int offset = -rhs[0].getIndex() + this->Expr[0].getIndex();
  rhsCopy.moveIndex(offset);
  return (*this) == rhsCopy;
}

//-------------------------------------------------FermiPolynomial<MonomialType>---------

template<typename MonomialType>
FermiPolynomial<MonomialType> & FermiPolynomial<MonomialType>::operator=(
    const FermiPolynomial<MonomialType> & rhs) {
  if (this != &rhs) {
    this->Terms = rhs.Terms;
  }
  return *this;
}

template<typename MonomialType>
bool FermiPolynomial<MonomialType>::isNorm() const {
  size_t len = this->Terms.size();
  for (size_t i = 0; i < len; i++) {
    if (!(this->Terms[i].second.isNorm())) {
      return false;
    }
  }
  return true;
}

template<typename MonomialType>
int FermiPolynomial<MonomialType>::findNonNorm() const {
  for (size_t i = 0; i < this->Terms.size(); i++) {
    if (!(this->Terms[i].second.isNorm())) {
      return i;
    }
  }
  return -1;
}

template<typename MonomialType>
void FermiPolynomial<MonomialType>::normOneTerm(int index) {
  //std::cout << "\nCurrent term: "
  //          << "index = " << index << ": " << this->Terms[index].second.toString()
  //          << std::endl;
  std::complex<double> pref = this->Terms[index].first;
  MonomialType mn(this->Terms[index].second);
  this->Terms.erase(this->Terms.begin() + index);
  (*this) += NormOnce(pref, mn);
}

template<typename MonomialType>
void FermiPolynomial<MonomialType>::normalize() {
  int index = -1;
  while ((index = findNonNorm()) >= 0) {
    normOneTerm(index);
  }
  this->eraseZeros();
}

template<typename MonomialType>
bool isNonNorm(std::pair<std::complex<double>, MonomialType> term) {
  return !(term.second.isNorm());
}

template<typename MonomialType>
void FermiPolynomial<MonomialType>::eraseNonNorm() {
  this->Terms.erase(
      std::remove_if(this->Terms.begin(), this->Terms.end(), isNonNorm<MonomialType>),
      this->Terms.end());
}

//---------------------------------------------------------Algebra Functions-------------

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > FermiCommute(const OpType & op1,
                                                     const OpType & op2) {
  FermiMonomial<OpType> mnr(op2);  //op2*op1
  mnr *= op1;
  if (op1.getCreatorF() == op2.getCreatorF()) {  //-op2*op1
    FermiPolynomial<FermiMonomial<OpType> > ans(std::complex<double>(-1, 0), mnr);
    return ans;
  }
  else {  //delta_{ij}-op2*op1
    FermiPolynomial<FermiMonomial<OpType> > ans(std::complex<double>(-1, 0), mnr);
    if (op1.getIndex() == op2.getIndex()) {
      OpType unit(true);
      FermiMonomial<OpType> mnu(unit);
      ans += mnu;
    }
    return ans;
  }
}

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > NormOnce(std::complex<double> pref,
                                                 FermiMonomial<OpType> mn) {
  size_t index = mn.findWrongOrder();
  FermiPolynomial<FermiMonomial<OpType> > mid = FermiCommute(mn[index], mn[index + 1]);
  FermiPolynomial<FermiMonomial<OpType> > ans;
  if (index > 0) {  //Wrong order not at beginning
    FermiMonomial<OpType> mnFront(mn.sliceExprStart(index));
    FermiPolynomial<FermiMonomial<OpType> > polyFront(std::complex<double>(1, 0),
                                                      mnFront);
    ans = mnFront;
    ans *= mid;
  }
  else {  //Wrong order at beginning
    ans = mid;
  }
  if (index < mn.getSize() - 2) {  //Not at end
    FermiPolynomial<FermiMonomial<OpType> > mnRear(std::complex<double>(1, 0),
                                                   mn.sliceExprEnd(index + 2));
    ans *= mnRear;
  }
  ans *= pref;
  return ans;
}

#endif  //QM_FERMI_OPERATORS_TEM_HPP
