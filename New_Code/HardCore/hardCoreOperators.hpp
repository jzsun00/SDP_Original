/*
  Jiazheng Sun
  Updated: May 14, 2024

  Class:
  HardCoreLadderOp, HardCore1DLadderOp, HardCoreMonomial, HardCorePolynomial.

  Function:

  Define ladder operators, monomials and polynomials for hard core Boson systems.
*/

#ifndef ORI_SDP_GS_HARDCOREOPERATORS_HPP
#define ORI_SDP_GS_HARDCOREOPERATORS_HPP

#include "../Basics/operators.hpp"
#include "../Basics/operators_Tem.cpp"

//---------------------------------------------------------------HardCoreLadderOp--------

template<typename IndexType>
class HardCoreLadderOp : public LadderOp<IndexType> {
 public:
  /*The constructors are identical to LadderOp.*/
  HardCoreLadderOp() : LadderOp<IndexType>() {}
  HardCoreLadderOp(IndexType index, bool creatorF) :
      LadderOp<IndexType>(index, creatorF) {}
  HardCoreLadderOp(bool isUnit) : LadderOp<IndexType>(isUnit) {}
  HardCoreLadderOp(HardCoreLadderOp const & rhs) : LadderOp<IndexType>(rhs) {}
  ~HardCoreLadderOp() {}
  /*Get information of the hard core Boson ladder operator.*/
  virtual std::string indexToString() const = 0;
  /*Overload operators.*/
  HardCoreLadderOp & operator=(HardCoreLadderOp const & rhs);
  bool operator==(HardCoreLadderOp const & rhs) const;
  bool operator!=(HardCoreLadderOp const & rhs) const { return !(*this == rhs); }
  virtual bool operator<(LadderOp<IndexType> const & rhs) const = 0;
  bool operator>(HardCoreLadderOp const & rhs) const {
    return !(*this < rhs || *this == rhs);
  }
};

//-------------------------------------------------------------HardCore1DLadderOp--------

class HardCore1DLadderOp : public HardCoreLadderOp<int> {
 public:
  /*The constructors are identical to HardCoreLadderOp.*/
  HardCore1DLadderOp() : HardCoreLadderOp<int>() {}
  HardCore1DLadderOp(int index, bool creatorF) : HardCoreLadderOp<int>(index, creatorF) {}
  HardCore1DLadderOp(bool isUnit) : HardCoreLadderOp<int>(isUnit) {}
  HardCore1DLadderOp(HardCore1DLadderOp const & rhs) : HardCoreLadderOp<int>(rhs) {}
  ~HardCore1DLadderOp() {}
  /*Get information of the 1-D hard core Boson ladder operator.*/
  virtual std::string indexToString() const { return std::to_string(this->index); }
  /*Overload operators.*/
  // If same type, compare index; if not same type, creator is always larger.
  virtual bool operator<(LadderOp<int> const & rhs) const;
  void moveIndex(int i) { this->index += i; }
};

//---------------------------------------------------------------HardCoreMonomial--------

template<typename OpType>
class HardCoreMonomial : public Monomial<OpType> {
 public:
  /*The constructors are identical to Monomial.*/
  HardCoreMonomial() : Monomial<OpType>() {}
  HardCoreMonomial(OpType & Op) : Monomial<OpType>(Op) {}
  HardCoreMonomial(vector<OpType> & Expr) : Monomial<OpType>(Expr) {}
  HardCoreMonomial(Monomial<OpType> const & rhs) : Monomial<OpType>(rhs) {}
  ~HardCoreMonomial() {}
  /*Tools for normalization.*/
  void reverse();
  int findWrongOrder() const;
  bool isNorm() const;
  HardCoreMonomial<OpType> sliceExprS(size_t index);
  HardCoreMonomial<OpType> sliceExprE(size_t index);
  /*Tools for infinite chain.*/
  void moveIndex(int i);
  bool equiv(HardCoreMonomial<OpType> const & rhs) const;
};

//--------------------------------------------------------------HardCorePolynomial-------

template<typename MonomialType>
class HardCorePolynomial : public Polynomial<MonomialType> {
 public:
  /*The constructors are identical to Polynomial.*/
  HardCorePolynomial() : Polynomial<MonomialType>() {}
  HardCorePolynomial(MonomialType const & mn) : Polynomial<MonomialType>(mn) {}
  HardCorePolynomial(complex<double> pref, MonomialType const & mn) :
      Polynomial<MonomialType>(pref, mn) {}
  HardCorePolynomial(HardCorePolynomial const & rhs) : Polynomial<MonomialType>(rhs) {}
  ~HardCorePolynomial() {}
  HardCorePolynomial & operator=(HardCorePolynomial<MonomialType> const & rhs);
  /*Tools for normalization.*/
  bool isNorm() const;
  int findNonNorm() const;
  void normOneTerm(int index);
  void normalize();
  void eraseNonNorm();
};

//----------------------------------------------------------------Algebra Functions------

template<typename OpType>
HardCorePolynomial<HardCoreMonomial<OpType> > HardCoreCommute(OpType op1, OpType op2);

template<typename OpType>
HardCorePolynomial<HardCoreMonomial<OpType> > HardCoreNormOnce(
    complex<double> pref,
    HardCoreMonomial<OpType> mn);

#include "./hardCoreOperators_Tem.cpp"

#endif  //ORI_SDP_GS_HARDCOREOPERATORS_HPP
