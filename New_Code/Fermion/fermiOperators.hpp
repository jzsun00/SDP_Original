/*
  Jiazheng Sun
  Updated: Jul 23, 2024

  Class:
  FermiLadderOp<IndexType>
  Fermi1DLadderOp
  FermiMonomial<OpType>
  FermiPolynomial<MonomialType>
  
  Define ladder operators, monomials and polynomials for Fermi systems.
*/

#ifndef QM_FERMI_OPERATORS_HPP
#define QM_FERMI_OPERATORS_HPP

#include "../Basics/operators.hpp"
#include "../Basics/operators_Tem.cpp"

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
  virtual bool operator<(LadderOp<IndexType> const & rhs) const = 0;
  bool operator>(FermiLadderOp const & rhs) const;
};

//------------------------------------------------------------Fermi1DLadderOp------------

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

//---------------------------------------------------------FermiMonomial<OpType>---------

template<typename OpType>
class FermiMonomial : public Monomial<OpType> {
 public:
  /*The constructors are identical to Monomial.*/
  FermiMonomial() : Monomial<OpType>() {}
  FermiMonomial(OpType & Op) : Monomial<OpType>(Op) {}
  FermiMonomial(std::vector<OpType> & Expr) : Monomial<OpType>(Expr) {}
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

//-------------------------------------------------FermiPolynomial<MonomialType>---------

template<typename MonomialType>
class FermiPolynomial : public Polynomial<MonomialType> {
 public:
  /*The constructors are identical to Polynomial.*/
  FermiPolynomial() : Polynomial<MonomialType>() {}
  FermiPolynomial(MonomialType const & mn) : Polynomial<MonomialType>(mn) {}
  FermiPolynomial(std::complex<double> pref, MonomialType const & mn) :
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

//---------------------------------------------------------Algebra Functions-------------

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > FermiCommute(OpType op1, OpType op2);

template<typename OpType>
FermiPolynomial<FermiMonomial<OpType> > NormOnce(std::complex<double> pref,
                                                 FermiMonomial<OpType> mn);

#endif  //QM_FERMI_OPERATORS_HPP
