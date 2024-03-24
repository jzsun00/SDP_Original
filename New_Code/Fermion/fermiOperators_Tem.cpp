/*
  Jiazheng Sun
  Updated: Mar 20, 2024

  Implementations of methods in class:
  Fermi1DLadderOp, FermiMonomial, FermiPolynomial.
 */

#ifndef ORI_SDP_GS_FERMIOPERATORS_TEM_CPP
#define ORI_SDP_GS_FERMIOPERATORS_TEM_CPP

#include "fermiOperators.hpp"

//------------------------------------------------------------------FermiLadderOp--------

template<typename IndexType>
FermiLadderOp<IndexType> & FermiLadderOp<IndexType>::operator=(
    FermiLadderOp<IndexType> const & rhs) {
  this->index = rhs.index;
  this->isUnit = rhs.isUnit;
  this->creatorF = rhs.creatorF;
  return *this;
}

template<typename IndexType>
bool FermiLadderOp<IndexType>::operator==(FermiLadderOp<IndexType> const & rhs) const {
  if (this->isUnit) {
    if (rhs.isUnit) {
      return true;
    }
    return false;
  }
  if (rhs.isUnit) {
    if (this->isUnit) {
      return true;
    }
    return false;
  }
  return (this->index == rhs.index) && (this->creatorF == rhs.creatorF);
}

//------------------------------------------------------------------FermiMonomial--------

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
FermiMonomial<OpType> FermiMonomial<OpType>::sliceExprS(size_t index) {
  return FermiMonomial<OpType>(
      vector<OpType>(this->Expr.begin(), this->Expr.begin() + index));
}

template<typename OpType>
FermiMonomial<OpType> FermiMonomial<OpType>::sliceExprE(size_t index) {
  return FermiMonomial<OpType>(
      vector<OpType>(this->Expr.begin() + index, this->Expr.end()));
}

//-----------------------------------------------------------------FermiPolynomial-------

template<typename MonomialType>
FermiPolynomial<MonomialType> & FermiPolynomial<MonomialType>::operator=(
    FermiPolynomial<MonomialType> const & rhs) {
  this->Terms = rhs.Terms;
  return *this;
}
/*
template<typename MonomialType>
void FermiPolynomial<MonomialType>::normalize() {
  for (typename vector<pair<complex<double>, MonomialType> >::iterator it =
           this->Terms.begin();
       it != this->Terms.end();
       ++it) {
    std::cout << "Length = " << this->Terms.size() << std::endl;
    std::cout << "it->second = " << it->second.toString() << std::endl;
    if (it->second.isNorm()) {
      continue;
    }
    std::cout << "NormOnce(it->first, it->second) = "
              << NormOnce(it->first, it->second).toString() << std::endl;
    (*this) += NormOnce(it->first, it->second);
    this->eraseZeros();
    std::cout << "(*this) = " << this->toString() << std::endl;
  }
  //this->eraseZeros();
}
*/

template<typename MonomialType>
bool FermiPolynomial<MonomialType>::isNorm() const {
  for (size_t i = 0; i < this->Terms.size(); i++) {
    if (!(this->Terms[i].second.isNorm())) {
      return false;
    }
  }
  return true;
}
/*
template<typename MonomialType>
void FermiPolynomial<MonomialType>::normalize() {
  for (size_t i = 0; i < this->Terms.size(); i++) {
    //std::cout << "Length = " << this->Terms.size() << std::endl;
    //std::cout << "\nCurrent term: "
    //          << "i = " << i << ": " << this->Terms[i].second.toString() << std::endl;
    if (this->Terms[i].second.isNorm()) {
      continue;
    }
    //std::cout << "NormOnce(it->first, it->second) = "
    //<< NormOnce(this->Terms[i].first, this->Terms[i].second).toString()
    //          << std::endl;
    (*this) += NormOnce(this->Terms[i].first, this->Terms[i].second);
    this->eraseZeros();
    std::cout << "\n(*this) = " << this->toString() << std::endl;
  }
  this->eraseZeros();
}
*/

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
  complex<double> pref = this->Terms[index].first;
  MonomialType mn(this->Terms[index].second);
  this->Terms.erase(this->Terms.begin() + index);
  (*this) += NormOnce(pref, mn);
  //std::cout << "\n(*this) = " << this->toString() << std::endl;
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
bool isNonNorm(pair<complex<double>, MonomialType> term) {
  return !(term.second.isNorm());
}

template<typename MonomialType>
void FermiPolynomial<MonomialType>::eraseNonNorm() {
  this->Terms.erase(
      std::remove_if(this->Terms.begin(), this->Terms.end(), isNonNorm<MonomialType>),
      this->Terms.end());
}

//----------------------------------------------------------------Algebra Functions------

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > FermiCommute(OpType op1, OpType op2) {
  FermiMonomial<OpType> mn0(op1);
  mn0 *= op2;
  FermiMonomial<OpType> mnr(op2);
  mnr *= op1;
  if (op1.getCreatorF() == op2.getCreatorF()) {
    FermiPolynomial<FermiMonomial<OpType> > ans(complex<double>(-1, 0), mnr);
    return ans;
  }
  else {
    FermiPolynomial<FermiMonomial<OpType> > ans(complex<double>(-1, 0), mnr);
    if (op1.getIndex() == op2.getIndex()) {
      OpType unit(true);
      FermiMonomial<OpType> mnu(unit);
      ans += mnu;
    }
    return ans;
  }
}

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > NormOnce(complex<double> pref,
                                                 FermiMonomial<OpType> mn) {
  size_t index = mn.findWrongOrder();
  //std::cout << "index = " << index << std::endl;
  FermiPolynomial<FermiMonomial<OpType> > mid = FermiCommute(mn[index], mn[index + 1]);
  //std::cout << "mid = " << mid.toString() << std::endl;
  FermiPolynomial<FermiMonomial<OpType> > ans;
  if (index > 0) {
    FermiMonomial<OpType> mnFront(mn.sliceExprS(index));
    FermiPolynomial<FermiMonomial<OpType> > polyFront(complex<double>(1, 0), mnFront);
    ans = mnFront;
    ans *= mid;
  }
  else {
    ans = mid;
  }
  //std::cout << "ans = " << ans.toString() << std::endl;
  if (index < mn.getSize() - 2) {
    FermiPolynomial<FermiMonomial<OpType> > mnRear(complex<double>(1, 0),
                                                   mn.sliceExprE(index + 2));
    //std::cout << "Rear = " << mnRear.toString() << std::endl;
    ans *= mnRear;
  }
  //std::cout << "ans = " << ans.toString() << std::endl;
  ans *= pref;
  //std::cout << "ans = " << ans.toString() << std::endl;
  return ans;
}

#endif  //ORI_SDP_GS_FERMIOPERATORS_TEM_CPP
