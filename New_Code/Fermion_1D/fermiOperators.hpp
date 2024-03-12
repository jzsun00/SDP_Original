/*
  Jiazheng Sun
  Updated: Mar 9, 2024

  Define ladder operators, monomials and polynomials for Fermi systems.
*/

#ifndef ORI_SDP_GS_FERMIOPERATORS_HPP
#define ORI_SDP_GS_FERMIOPERATORS_HPP

#include "./fermiStates.hpp"
#include "./operators.hpp"

//------------------------------------------------------------------FermiLadderOp--------

class FermiLadderOp : public LadderOp {
 public:
  /*The constructors are identical to LadderOp.*/
  FermiLadderOp(int index, bool creatorF) : LadderOp(index, creatorF) {}
  FermiLadderOp(LadderOp const & rhs) : LadderOp(rhs) {}
  /*Define operators at Fock states.*/
  FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
};

//------------------------------------------------------------------FermiMonomial--------

class FermiMonomial : public Monomial<FermiLadderOp> {
 public:
  /*The constructors are identical to Monomial.*/
  FermiMonomial() : Monomial<FermiLadderOp>() {}
  FermiMonomial(FermiLadderOp & Op) : Monomial<FermiLadderOp>(Op) {}
  FermiMonomial(Monomial<FermiLadderOp> const & rhs) : Monomial<FermiLadderOp>(rhs) {}
  ~FermiMonomial() {}
  /*Define operators at Fock states.*/
  FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
};

//-----------------------------------------------------------------FermiPolynomial-------

class FermiPolynomial : public Polynomial<FermiLadderOp> {
 public:
  /*The constructors are identical to Polynomial.*/
  FermiPolynomial() : Polynomial<FermiLadderOp>() {}
  FermiPolynomial(Monomial<FermiLadderOp> const & mn) : Polynomial<FermiLadderOp>(mn) {}
  FermiPolynomial(FermiPolynomial const & rhs) : Polynomial<FermiLadderOp>(rhs) {}
  ~FermiPolynomial() {}
  /*Define operators at Fock states.*/
  FermiState operator*(FermiFockState const & rhs) const;
  //FermiState operator*(FermiState const & rhs) const;
};

#endif  //ORI_SDP_GS_FERMIOPERATORS_HPP
