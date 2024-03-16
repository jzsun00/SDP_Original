/*
  Jiazheng Sun
  Updated: Mar 16, 2024

  Class:
  SpinHalfOp, SpinHalfMonomial, SpinHalfPolynomial.

  Define operators, monomials and polynomials for spin-1/2 systems.
*/

#ifndef ORI_SDP_GS_SPINOPERATORS1D_HPP
#define ORI_SDP_GS_SPINOPERATORS1D_HPP

#include "../Basics/operators.hpp"
#include "./spinStates1D.hpp"

//--------------------------------------------------------------- ---SpinHalfOp----------

class SpinHalfOp : public SpinOp {
 public:
  /*Construct a spin-1/2 operator with specified index and type,
   default constructor use INT_MIN and z-component.*/
  SpinHalfOp() : SpinOp() {}
  SpinHalfOp(int index) : SpinOp(index) {}
  SpinHalfOp(int index, bool isPlus) : SpinOp(index, isPlus) {}
  SpinHalfOp(int index, bool isZ, bool isPlus) : SpinOp(index, isZ, isPlus) {}
  SpinHalfOp(SpinHalfOp const & rhs) : SpinOp(rhs) {}
  ~SpinHalfOp() {}
  /*Overload operators.*/
  SpinHalfState operator*(SpinHalfBaseState const & rhs) const;
  SpinHalfState operator*(SpinHalfState const & rhs) const;
};

//----------------------------------------------------------------SpinHalfMonomial-------

class SpinHalfMonomial : public Monomial<SpinHalfOp> {
 public:
  SpinHalfMonomial() : Monomial<SpinHalfOp>() {}
  SpinHalfMonomial(SpinHalfOp & Op) : Monomial<SpinHalfOp>(Op) {}
  SpinHalfMonomial(SpinHalfMonomial const & rhs) : Monomial<SpinHalfOp>(rhs) {}
  ~SpinHalfMonomial() {}
  /*Overload operators.*/
  SpinHalfState operator*(SpinHalfBaseState const & rhs) const;
  SpinHalfState operator*(SpinHalfState const & rhs) const;
};

//---------------------------------------------------------------SpinHalfPolynomial------

class SpinHalfPolynomial : public Polynomial<SpinHalfMonomial> {
 public:
  SpinHalfPolynomial() : Polynomial<SpinHalfMonomial>() {}
  SpinHalfPolynomial(SpinHalfMonomial const & mn) : Polynomial<SpinHalfMonomial>(mn) {}
  SpinHalfPolynomial(SpinHalfPolynomial const & rhs) :
      Polynomial<SpinHalfMonomial>(rhs) {}
  ~SpinHalfPolynomial() {}
  /*Overload operators.*/
  SpinHalfState operator*(SpinHalfBaseState const & rhs) const;
  SpinHalfState operator*(SpinHalfState const & rhs) const;
};

#endif  //ORI_SDP_GS_SPINOPERATORS1D_HPP
