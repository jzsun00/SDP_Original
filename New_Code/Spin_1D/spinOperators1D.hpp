/*
  Jiazheng Sun
  Updated: Aug 2, 2024

  Class:
  SpinHalfOp1D
  SpinHalfMonomial1D
  SpinHalfPolynomial1D

  Define operators, monomials and polynomials for 1D spin systems.
*/

#ifndef QM_SPIN_OPERATORS_1D_HPP
#define QM_SPIN_OPERATORS_1D_HPP

#include "../Basics/operators_Tem.hpp"
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
  SpinHalfOp1D(const SpinHalfOp1D & rhs) : SpinOp<int>(rhs) {}
  virtual ~SpinHalfOp1D() {}
  /*Get information of the spin operator.*/
  virtual std::string indexToString() const { return std::to_string(index); }
  /*Overload operators.*/
  SpinHalfOp1D & operator=(const SpinHalfOp1D & rhs);
  SpinHalfState1D operator*(const SpinHalfBaseState1D & rhs) const;
  SpinHalfState1D operator*(const SpinHalfState1D & rhs) const;
};

//-----------------------------------------------------------SpinHalfMonomial1D----------

class SpinHalfMonomial1D : public Monomial<SpinHalfOp1D> {
 public:
  SpinHalfMonomial1D() : Monomial<SpinHalfOp1D>() {}
  SpinHalfMonomial1D(SpinHalfOp1D & Op) : Monomial<SpinHalfOp1D>(Op) {}
  SpinHalfMonomial1D(SpinHalfMonomial1D const & rhs) : Monomial<SpinHalfOp1D>(rhs) {}
  ~SpinHalfMonomial1D() {}
  /*Overload operators.*/
  SpinHalfMonomial1D & operator=(const SpinHalfMonomial1D & rhs);
  SpinHalfState1D operator*(const SpinHalfBaseState1D & rhs) const;
  SpinHalfState1D operator*(const SpinHalfState1D & rhs) const;
};

//----------------------------------------------------------SpinHalfPolynomial1D---------

class SpinHalfPolynomial1D : public Polynomial<SpinHalfMonomial1D> {
 public:
  SpinHalfPolynomial1D() : Polynomial<SpinHalfMonomial1D>() {}
  SpinHalfPolynomial1D(const SpinHalfMonomial1D & mn) :
      Polynomial<SpinHalfMonomial1D>(mn) {}
  SpinHalfPolynomial1D(const SpinHalfPolynomial1D & rhs) :
      Polynomial<SpinHalfMonomial1D>(rhs) {}
  ~SpinHalfPolynomial1D() {}
  /*Overload operators.*/
  SpinHalfPolynomial1D & operator=(const SpinHalfPolynomial1D & rhs);
  SpinHalfState1D operator*(const SpinHalfBaseState1D & rhs) const;
  SpinHalfState1D operator*(const SpinHalfState1D & rhs) const;
  void operate(const SpinHalfState1D & rhs, SpinHalfState1D & result) const;
};

#endif  //QM_SPIN_OPERATORS_1D_HPP
