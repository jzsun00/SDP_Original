/*
  Jiazheng Sun
  Updated: Jun 19, 2024

  Class:
  SpinHalfOp1D
  SpinHalfMonomial1D
  SpinHalfPolynomial1D

  Define operators, monomials and polynomials for 1D spin systems.
*/

#ifndef QM_SPINOPERATORS1D_HPP
#define QM_SPINOPERATORS1D_HPP

#include "../Basics/operators.hpp"
#include "../Basics/operators_Tem.cpp"
#include "./spinStates1D.hpp"

//--------------------------------------------------------------SpinHalfOp1D-------------

class SpinHalfOp1D : public SpinOp<int> {
 public:
  /*Construct a spin-1/2 operator with specified index and type,
   default constructor use INT_MIN and z-component.*/
  SpinHalfOp1D() : SpinOp<int>() {}
  SpinHalfOp1D(int index) : SpinOp<int>(index) {}
  SpinHalfOp1D(int index, bool isPlus) : SpinOp<int>(index, isPlus) {}
  SpinHalfOp1D(int index, bool isZ, bool isPlus) : SpinOp<int>(index, isZ, isPlus) {}
  SpinHalfOp1D(SpinHalfOp1D const & rhs) : SpinOp<int>(rhs) {}
  ~SpinHalfOp1D() {}
  /*Get information of the spin operator.*/
  virtual std::string indexToString() const { return std::to_string(index); }
  /*Overload operators.*/
  SpinHalfState1D operator*(SpinHalfBaseState1D const & rhs) const;
  SpinHalfState1D operator*(SpinHalfState1D const & rhs) const;
};

//-----------------------------------------------------------SpinHalfMonomial1D----------

class SpinHalfMonomial1D : public Monomial<SpinHalfOp1D> {
 public:
  SpinHalfMonomial1D() : Monomial<SpinHalfOp1D>() {}
  SpinHalfMonomial1D(SpinHalfOp1D & Op) : Monomial<SpinHalfOp1D>(Op) {}
  SpinHalfMonomial1D(SpinHalfMonomial1D const & rhs) : Monomial<SpinHalfOp1D>(rhs) {}
  ~SpinHalfMonomial1D() {}
  /*Overload operators.*/
  SpinHalfState1D operator*(SpinHalfBaseState1D const & rhs) const;
  SpinHalfState1D operator*(SpinHalfState1D const & rhs) const;
};

//----------------------------------------------------------SpinHalfPolynomial1D---------

class SpinHalfPolynomial1D : public Polynomial<SpinHalfMonomial1D> {
 public:
  SpinHalfPolynomial1D() : Polynomial<SpinHalfMonomial1D>() {}
  SpinHalfPolynomial1D(SpinHalfMonomial1D const & mn) :
      Polynomial<SpinHalfMonomial1D>(mn) {}
  SpinHalfPolynomial1D(Polynomial<SpinHalfMonomial1D> const & rhs) :
      Polynomial<SpinHalfMonomial1D>(rhs) {}
  ~SpinHalfPolynomial1D() {}
  /*Overload operators.*/
  SpinHalfPolynomial1D & operator=(SpinHalfPolynomial1D const & rhs);
  SpinHalfState1D operator*(SpinHalfBaseState1D const & rhs) const;
  SpinHalfState1D operator*(SpinHalfState1D const & rhs) const;
};

#endif  //QM_SPINOPERATORS1D_HPP
