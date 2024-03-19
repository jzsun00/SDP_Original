/*
  Jiazheng Sun
  Updated: Mar 18, 2024

  Class:
  Fermi1DLadderOp, FermiMonomial, FermiPolynomial

  Define ladder operators, monomials and polynomials for Fermi systems.
*/

#ifndef ORI_SDP_GS_FERMIOPERATORS_HPP
#define ORI_SDP_GS_FERMIOPERATORS_HPP

#include "../Basics/operators.hpp"
#include "./fermiStates.hpp"

//-----------------------------------------------------------------Fermi1DLadderOp-------

class Fermi1DLadderOp : public LadderOp<int> {
 public:
  /*The constructors are identical to LadderOp.*/
  Fermi1DLadderOp() : LadderOp<int>() {}
  Fermi1DLadderOp(int index, bool creatorF) : LadderOp<int>(index, creatorF) {}
  Fermi1DLadderOp(Fermi1DLadderOp const & rhs) : LadderOp<int>(rhs) {}
  ~Fermi1DLadderOp() {}
  /*Get information of the 1-D Fermi ladder operator.*/
  virtual std::string indexToString() const { return std::to_string(this->index); }
  /*Overload operators.*/
  // If same type, compare index; if not same type, creator is always larger.
  virtual bool operator<(LadderOp const & rhs) const;
  /*Define operators at Fock states.*/
  FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
};

//------------------------------------------------------------------FermiMonomial--------

template<typename OpType>
class FermiMonomial : public Monomial<OpType> {
 public:
  /*The constructors are identical to Monomial.*/
  FermiMonomial() : Monomial<OpType>() {}
  FermiMonomial(OpType & Op) : Monomial<OpType>(Op) {}
  FermiMonomial(Monomial<OpType> const & rhs) : Monomial<OpType>(rhs) {}
  ~FermiMonomial() {}
  /*Define operators at Fock states.*/
  FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
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
  /*Define operators at Fock states.*/
  FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
};

#endif  //ORI_SDP_GS_FERMIOPERATORS_HPP
