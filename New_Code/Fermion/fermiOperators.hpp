/*
  Jiazheng Sun
  Updated: Mar 21, 2024

  Class:
  FermiLadderOp, Fermi1DLadderOp, FermiMonomial, FermiPolynomial.

  Define ladder operators, monomials and polynomials for Fermi systems.
*/

#ifndef ORI_SDP_GS_FERMIOPERATORS_HPP
#define ORI_SDP_GS_FERMIOPERATORS_HPP

#include "../Basics/operators.hpp"
//#include "./fermiStates.hpp"

//------------------------------------------------------------------FermiLadderOp--------

template<typename IndexType>
class FermiLadderOp : public LadderOp<IndexType> {
 public:
  /*The constructors are identical to LadderOp.*/
  FermiLadderOp() : LadderOp<IndexType>() {}
  FermiLadderOp(IndexType index, bool creatorF) : LadderOp<IndexType>(index, creatorF) {}
  FermiLadderOp(bool isUnit) : LadderOp<IndexType>(isUnit) {}
  FermiLadderOp(FermiLadderOp const & rhs) : LadderOp<IndexType>(rhs) {}
  ~FermiLadderOp() {}
  /*Get information of the Fermi ladder operator.*/
  virtual std::string indexToString() const = 0;
  /*Overload operators.*/
  FermiLadderOp & operator=(FermiLadderOp const & rhs);
  bool operator==(FermiLadderOp const & rhs) const;
  bool operator!=(FermiLadderOp const & rhs) const { return !(*this == rhs); }
  virtual bool operator<(LadderOp<IndexType> const & rhs) const = 0;
  bool operator>(FermiLadderOp const & rhs) const {
    return !(*this < rhs || *this == rhs);
  }
};

//-----------------------------------------------------------------Fermi1DLadderOp-------

class Fermi1DLadderOp : public FermiLadderOp<int> {
 public:
  /*The constructors are identical to FermiLadderOp.*/
  Fermi1DLadderOp() : FermiLadderOp<int>() {}
  Fermi1DLadderOp(int index, bool creatorF) : FermiLadderOp<int>(index, creatorF) {}
  Fermi1DLadderOp(bool isUnit) : FermiLadderOp<int>(isUnit) {}
  Fermi1DLadderOp(Fermi1DLadderOp const & rhs) : FermiLadderOp<int>(rhs) {}
  ~Fermi1DLadderOp() {}
  /*Get information of the 1-D Fermi ladder operator.*/
  virtual std::string indexToString() const { return std::to_string(this->index); }
  /*Overload operators.*/
  // If same type, compare index; if not same type, creator is always larger.
  virtual bool operator<(LadderOp<int> const & rhs) const;
  /*Define operators at Fock states.*/
  //FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
};

//------------------------------------------------------------------FermiMonomial--------

template<typename OpType>
class FermiMonomial : public Monomial<OpType> {
 public:
  /*The constructors are identical to Monomial.*/
  FermiMonomial() : Monomial<OpType>() {}
  FermiMonomial(OpType & Op) : Monomial<OpType>(Op) {}
  FermiMonomial(vector<OpType> & Expr) : Monomial<OpType>(Expr) {}
  FermiMonomial(Monomial<OpType> const & rhs) : Monomial<OpType>(rhs) {}
  ~FermiMonomial() {}
  /*Define operators at Fock states.*/
  //FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
  /*Tools for normalization.*/
  int findWrongOrder() const;
  bool isNorm() const;
  FermiMonomial<OpType> sliceExprS(size_t index);
  FermiMonomial<OpType> sliceExprE(size_t index);
};

//-----------------------------------------------------------------FermiPolynomial-------

template<typename MonomialType>
class FermiPolynomial : public Polynomial<MonomialType> {
 public:
  /*The constructors are identical to Polynomial.*/
  FermiPolynomial() : Polynomial<MonomialType>() {}
  FermiPolynomial(MonomialType const & mn) : Polynomial<MonomialType>(mn) {}
  FermiPolynomial(complex<double> pref, MonomialType const & mn) :
      Polynomial<MonomialType>(pref, mn) {}
  FermiPolynomial(FermiPolynomial const & rhs) : Polynomial<MonomialType>(rhs) {}
  ~FermiPolynomial() {}
  FermiPolynomial & operator=(FermiPolynomial<MonomialType> const & rhs);
  /*Define operators at Fock states.*/
  //FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
  /*Tools for normalization.*/
  bool isNorm() const;
  int findNonNorm() const;
  void normOneTerm(int index);
  void normalize();
  void eraseNonNorm();
};

//----------------------------------------------------------------Algebra Functions------

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > FermiCommute(OpType op1, OpType op2);

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > NormOnce(complex<double> pref,
                                                 FermiMonomial<OpType> mn);

#include "./fermiOperators_Tem.cpp"

#endif  //ORI_SDP_GS_FERMIOPERATORS_HPP
