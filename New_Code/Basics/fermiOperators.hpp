/*
  Jan 31, 2024  Jiazheng Sun

  Define ladder operators, monomials and polynomials.
*/

#ifndef ORI_SDP_GS_FERMIOPERATORS_HPP
#define ORI_SDP_GS_FERMIOPERATORS_HPP

#include "./fermiStates.hpp"
#include "./operators.hpp"

class FermiLadderOp : public LadderOp {
 public:
  FermiLadderOp() : LadderOp() {}
  FermiLadderOp(int index, bool creatorF) : LadderOp(index, creatorF) {}
  FermiLadderOp(LadderOp const & rhs) : LadderOp(rhs) {}
  FermiState operator*(FermiFockState const & rhs) const;
  FermiState operator*(FermiState const & rhs) const;
};

class FermiMonomial : public Monomial {
 public:
  FermiMonomial(Monomial const & rhs) : Monomial(rhs) {}
  FermiState operator*(FermiFockState const & rhs) const;
  FermiState operator*(FermiState const & rhs) const;
};

class FermiPolynomial : public Polynomial {
 public:
  FermiState operator*(FermiFockState const & rhs) const;
  FermiState operator*(FermiState const & rhs) const;
};

#endif
