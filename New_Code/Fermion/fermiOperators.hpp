/*
  Jiazheng Sun
  Updated: Mar 18, 2024

  Class:
  Fermi1DLadderOp, FermiMonomial, FermiPolynomial

  Define ladder operators, monomials and polynomials for Fermi systems.
*/

#ifndef ORI_SDP_GS_FERMIOPERATORS_HPP
#define ORI_SDP_GS_FERMIOPERATORS_HPP

//#include "./fermiStates.hpp"
#include "../Basics/operators.hpp"

//-----------------------------------------------------------------Fermi1DLadderOp-------

class Fermi1DLadderOp : public LadderOp<int> {
 public:
  /*The constructors are identical to LadderOp.*/
  Fermi1DLadderOp() : LadderOp<int>() {}
  Fermi1DLadderOp(int index, bool creatorF) : LadderOp<int>(index, creatorF) {}
  Fermi1DLadderOp(Fermi1DLadderOp const & rhs) : LadderOp<int>(rhs) {}
  ~Fermi1DLadderOp() {}
  /*Get information of the 1-D Fermi ladder operator.*/
  virtual int getMinValue() const { return std::numeric_limits<int>::min(); }
  virtual std::string indexToString() const { return std::to_string(this->index); }
  /*Overload operators.*/
  virtual bool operator<(LadderOp const & rhs) const;
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
  FermiMonomial(Fermi1DLadderOp & Op) : Monomial<OpType>(Op) {}
  FermiMonomial(Monomial<OpType> const & rhs) : Monomial<OpType>(rhs) {}
  ~FermiMonomial() {}
  /*Define operators at Fock states.*/
  //FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
};

//-----------------------------------------------------------------FermiPolynomial-------

template<typename MonomialType>
class FermiPolynomial : public Polynomial<MonomialType> {
 public:
  /*The constructors are identical to Polynomial.*/
  FermiPolynomial() : Polynomial<MonomialType>() {}
  FermiPolynomial(Monomial<MonomialType> const & mn) : Polynomial<MonomialType>(mn) {}
  FermiPolynomial(FermiPolynomial const & rhs) : Polynomial<MonomialType>(rhs) {}
  ~FermiPolynomial() {}
  /*Define operators at Fock states.*/
  //FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
};

#endif  //ORI_SDP_GS_FERMIOPERATORS_HPP
