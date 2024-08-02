/*
  Jiazheng Sun
  Updated: Aug 2, 2024

  Class:
  FermiLadderOp<IndexType>
  Fermi1DLadderOp
  FermiMonomial<OpType>
  FermiPolynomial<MonomialType>
  
  Define ladder operators, monomials and polynomials for Fermi systems.
*/

#ifndef QM_FERMI_OPERATORS_HPP
#define QM_FERMI_OPERATORS_HPP

#include "../Basics/operators_Tem.hpp"

//------------------------------------------------------FermiLadderOp<IndexType>---------

template<typename IndexType>
class FermiLadderOp : public LadderOp<IndexType> {
 public:
  /*The constructors are identical to LadderOp.*/
  FermiLadderOp() : LadderOp<IndexType>() {}
  FermiLadderOp(IndexType index, bool creatorF) : LadderOp<IndexType>(index, creatorF) {}
  FermiLadderOp(bool isUnit) : LadderOp<IndexType>(isUnit) {}
  FermiLadderOp(const FermiLadderOp<IndexType> & rhs) : LadderOp<IndexType>(rhs) {}
  virtual ~FermiLadderOp() {}
  /*Get information of the Fermi ladder operator.*/
  virtual std::string indexToString() const = 0;
  /*Overload operators.*/
  FermiLadderOp<IndexType> & operator=(const FermiLadderOp<IndexType> & rhs);
  virtual bool operator<(const FermiLadderOp<IndexType> & rhs) const = 0;
  bool operator>(const FermiLadderOp & rhs) const;
};

//------------------------------------------------------------Fermi1DLadderOp------------

class Fermi1DLadderOp : public FermiLadderOp<int> {
 public:
  /*The constructors are identical to FermiLadderOp.*/
  Fermi1DLadderOp() : FermiLadderOp<int>() {}
  Fermi1DLadderOp(int index, bool creatorF) : FermiLadderOp<int>(index, creatorF) {}
  Fermi1DLadderOp(bool isUnit) : FermiLadderOp<int>(isUnit) {}
  Fermi1DLadderOp(const Fermi1DLadderOp & rhs) : FermiLadderOp<int>(rhs) {}
  virtual ~Fermi1DLadderOp() {}
  /*Get information of the 1-D Fermi ladder operator.*/
  virtual std::string indexToString() const { return std::to_string(this->index); }
  /*Overload operators.*/
  Fermi1DLadderOp & operator=(const Fermi1DLadderOp & rhs);
  /*If same type, compare index:
    annihilation operators ascending, creation operators descending;
    if not same type, creator is always larger.*/
  virtual bool operator<(const FermiLadderOp<int> & rhs) const;
  void moveIndex(int i) { this->index += i; }
};

//---------------------------------------------------------FermiMonomial<OpType>---------

template<typename OpType>
class FermiMonomial : public Monomial<OpType> {
 public:
  /*The constructors are identical to Monomial.*/
  FermiMonomial() : Monomial<OpType>() {}
  FermiMonomial(const OpType & Op) : Monomial<OpType>(Op) {}
  FermiMonomial(const std::vector<OpType> & Expr) : Monomial<OpType>(Expr) {}
  FermiMonomial(const FermiMonomial<OpType> & rhs) : Monomial<OpType>(rhs) {}
  virtual ~FermiMonomial() {}
  /*Overload operators.*/
  FermiMonomial<OpType> & operator=(const FermiMonomial<OpType> & rhs);
  /*Tools for normal order.
    Normal order: all creation operators on the right,
    annihilation operators index ascending, creation operators index descending.*/
  void reverse();
  int findWrongOrder() const;
  bool isNorm() const;
  FermiMonomial<OpType> sliceExprStart(size_t index) const;
  FermiMonomial<OpType> sliceExprEnd(size_t index) const;
  /*Tools for infinite system.*/
  void moveIndex(int i);
  bool equiv(const FermiMonomial<OpType> & rhs) const;
};

//-------------------------------------------------FermiPolynomial<MonomialType>---------

template<typename MonomialType>
class FermiPolynomial : public Polynomial<MonomialType> {
 public:
  /*The constructors are identical to Polynomial.*/
  FermiPolynomial() : Polynomial<MonomialType>() {}
  FermiPolynomial(const MonomialType & mn) : Polynomial<MonomialType>(mn) {}
  FermiPolynomial(std::complex<double> pref, const MonomialType & mn) :
      Polynomial<MonomialType>(pref, mn) {}
  FermiPolynomial(const FermiPolynomial & rhs) : Polynomial<MonomialType>(rhs) {}
  virtual ~FermiPolynomial() {}
  FermiPolynomial & operator=(const FermiPolynomial<MonomialType> & rhs);
  /*Tools for normalization.*/
  bool isNorm() const;
  int findNonNorm() const;
  void normOneTerm(int index);
  void normalize();
  void eraseNonNorm();
};

//---------------------------------------------------------Algebra Functions-------------

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > FermiCommute(const OpType & op1,
                                                     const OpType & op2);

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > NormOnce(std::complex<double> pref,
                                                 FermiMonomial<OpType> mn);

#endif  //QM_FERMI_OPERATORS_HPP
