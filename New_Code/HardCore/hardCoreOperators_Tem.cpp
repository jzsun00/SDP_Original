/*
  Jiazheng Sun
  Updated: Mar 20, 2024

  Implementations of methods in class:
  Fermi1DLadderOp, FermiMonomial, FermiPolynomial.
 */

#ifndef ORI_SDP_GS_HARDCOREOPERATORS_TEM_CPP
#define ORI_SDP_GS_HARDCOREOPERATORS_TEM_CPP

#include "hardCoreOperators.hpp"

//---------------------------------------------------------------HardCoreLadderOp--------

template<typename IndexType>
HardCoreLadderOp<IndexType> & HardCoreLadderOp<IndexType>::operator=(
    HardCoreLadderOp<IndexType> const & rhs) {
  this->index = rhs.index;
  this->isUnit = rhs.isUnit;
  this->creatorF = rhs.creatorF;
  return *this;
}

template<typename IndexType>
bool HardCoreLadderOp<IndexType>::operator==(
    HardCoreLadderOp<IndexType> const & rhs) const {
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

//---------------------------------------------------------------HardCoreMonomial--------

template<typename OpType>
int HardCoreMonomial<OpType>::findWrongOrder() const {
  for (size_t i = 0; i < this->Expr.size() - 1; i++) {
    if (this->Expr[i] > this->Expr[i + 1]) {
      return i;
    }
  }
  return -1;
}

template<typename OpType>
bool HardCoreMonomial<OpType>::isNorm() const {
  if (this->Expr.size() <= 1) {
    return true;
  }
  if (findWrongOrder() < 0) {
    return true;
  }
  return false;
}

template<typename OpType>
HardCoreMonomial<OpType> HardCoreMonomial<OpType>::sliceExprS(size_t index) {
  return HardCoreMonomial<OpType>(
      vector<OpType>(this->Expr.begin(), this->Expr.begin() + index));
}

template<typename OpType>
HardCoreMonomial<OpType> HardCoreMonomial<OpType>::sliceExprE(size_t index) {
  return HardCoreMonomial<OpType>(
      vector<OpType>(this->Expr.begin() + index, this->Expr.end()));
}

//--------------------------------------------------------------HardCorePolynomial-------

template<typename MonomialType>
HardCorePolynomial<MonomialType> & HardCorePolynomial<MonomialType>::operator=(
    HardCorePolynomial<MonomialType> const & rhs) {
  this->Terms = rhs.Terms;
  return *this;
}

template<typename MonomialType>
bool HardCorePolynomial<MonomialType>::isNorm() const {
  for (size_t i = 0; i < this->Terms.size(); i++) {
    if (!(this->Terms[i].second.isNorm())) {
      return false;
    }
  }
  return true;
}

template<typename MonomialType>
int HardCorePolynomial<MonomialType>::findNonNorm() const {
  for (size_t i = 0; i < this->Terms.size(); i++) {
    if (!(this->Terms[i].second.isNorm())) {
      return i;
    }
  }
  return -1;
}

template<typename MonomialType>
void HardCorePolynomial<MonomialType>::normOneTerm(int index) {
  //std::cout << "\nCurrent term: "
  //          << "index = " << index << ": " << this->Terms[index].second.toString()
  //          << std::endl;
  complex<double> pref = this->Terms[index].first;
  MonomialType mn(this->Terms[index].second);
  this->Terms.erase(this->Terms.begin() + index);
  (*this) += HardCoreNormOnce(pref, mn);
  //std::cout << "\n(*this) = " << this->toString() << std::endl;
}

template<typename MonomialType>
void HardCorePolynomial<MonomialType>::normalize() {
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
void HardCorePolynomial<MonomialType>::eraseNonNorm() {
  this->Terms.erase(
      std::remove_if(this->Terms.begin(), this->Terms.end(), isNonNorm<MonomialType>),
      this->Terms.end());
}

//----------------------------------------------------------------Algebra Functions------

template<typename OpType>
HardCorePolynomial<HardCoreMonomial<OpType> > HardCoreCommute(OpType op1, OpType op2) {
  HardCoreMonomial<OpType> mn0(op1);
  mn0 *= op2;
  HardCoreMonomial<OpType> mnr(op2);
  mnr *= op1;
  if (op1.getCreatorF() == op2.getCreatorF()) {
    //HardCorePolynomial<HardCoreMonomial<OpType> > ans(complex<double>(1, 0), mnr);
    HardCorePolynomial<HardCoreMonomial<OpType> > ans(complex<double>(-1, 0), mnr);
    return ans;
  }
  if (op1.getIndex() != op2.getIndex()) {
    //HardCorePolynomial<HardCoreMonomial<OpType> > ans(complex<double>(1, 0), mnr);
    HardCorePolynomial<HardCoreMonomial<OpType> > ans(complex<double>(-1, 0), mnr);
    return ans;
  }
  else {
    HardCorePolynomial<HardCoreMonomial<OpType> > ans(complex<double>(-1, 0), mnr);
    OpType unit(true);
    HardCoreMonomial<OpType> mnu(unit);
    ans += mnu;
    return ans;
  }
}

template<typename OpType>
HardCorePolynomial<HardCoreMonomial<OpType> > HardCoreNormOnce(
    complex<double> pref,
    HardCoreMonomial<OpType> mn) {
  size_t index = mn.findWrongOrder();
  //std::cout << "index = " << index << std::endl;
  HardCorePolynomial<HardCoreMonomial<OpType> > mid =
      HardCoreCommute(mn[index], mn[index + 1]);
  //std::cout << "mid = " << mid.toString() << std::endl;
  HardCorePolynomial<HardCoreMonomial<OpType> > ans;
  if (index > 0) {
    HardCoreMonomial<OpType> mnFront(mn.sliceExprS(index));
    HardCorePolynomial<HardCoreMonomial<OpType> > polyFront(complex<double>(1, 0),
                                                            mnFront);
    ans = mnFront;
    ans *= mid;
  }
  else {
    ans = mid;
  }
  //std::cout << "ans = " << ans.toString() << std::endl;
  if (index < mn.getSize() - 2) {
    HardCorePolynomial<HardCoreMonomial<OpType> > mnRear(complex<double>(1, 0),
                                                         mn.sliceExprE(index + 2));
    //std::cout << "Rear = " << mnRear.toString() << std::endl;
    ans *= mnRear;
  }
  //std::cout << "ans = " << ans.toString() << std::endl;
  ans *= pref;
  //std::cout << "ans = " << ans.toString() << std::endl;
  return ans;
}

#endif  //ORI_SDP_GS_HARDCOREOPERATORS_TEM_CPP
