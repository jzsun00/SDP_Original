/*
  Jiazheng Sun
  Updated: Mar 12, 2024

  Define spin operators, monomials and polynomials for spin systems.
*/

#ifndef ORI_SDP_GS_SPINOPERATORS1D_HPP
#define ORI_SDP_GS_SPINOPERATORS1D_HPP

#include "../Basics/operators.hpp"
#include "./mixMonomial.hpp"
#include "./spinStates1D.hpp"

#define PLUS 0
#define MINUS 1
#define SZ 2

//--------------------------------------------------------------------SpinHalfOp---------

class SpinHalfOp {
 protected:
  int index;

 public:
  /*Construct a spin operator with specified index and type,
   default constructor use INT_MIN and z-component.*/
  SpinHalfOp() : index(std::numeric_limits<int>::min()) {}
  SpinHalfOp(int index) : index(index) {}
  SpinHalfOp(SpinHalfOp const & rhs) { index = rhs.index; }
  virtual ~SpinHalfOp() {}
  /*Get information of the spin operator.*/
  int getIndex() const { return index; }
  virtual std::string toString() const = 0;
  /*Overload operators.*/
  //SpinHalfOp & operator=(SpinOp const & rhs) {
  //  index = rhs.index;
  //  return *this;
  //}
  //bool operator==(SpinOp const & rhs) const { return index == rhs.index; };
  // Throw std::invalid_argument exception if operators are not the same kind
  //bool operator<(SpinOp const & rhs) const;
  //bool operator>(SpinOp const & rhs) const { return !(*this < rhs || *this == rhs); }
  virtual void herm() = 0;
  virtual SpinHalfState operator*(SpinHalfBaseState const & rhs) const = 0;
  virtual SpinHalfState operator*(SpinHalfState const & rhs) const = 0;
};

//-------------------------------------------------------------------SpinZHalfOp---------

class SpinZHalfOp : public SpinHalfOp {
 public:
  SpinZHalfOp() : SpinHalfOp() {}
  SpinZHalfOp(int index) : SpinHalfOp(index) {}
  SpinZHalfOp(SpinHalfOp const & rhs) : SpinHalfOp(rhs) {}
  virtual ~SpinZHalfOp() {}
  virtual std::string toString() const;
  virtual void herm(){};
  virtual SpinHalfState operator*(SpinHalfBaseState const & rhs) const;
  virtual SpinHalfState operator*(SpinHalfState const & rhs) const;
};

//-------------------------------------------------------------------SpinUDHalfOp--------

class SpinUDHalfOp : public SpinHalfOp {
 private:
  bool plusF;

 public:
  SpinUDHalfOp() : SpinHalfOp(), plusF(false) {}
  SpinUDHalfOp(int index, bool plusF) : SpinHalfOp(index), plusF(plusF) {}
  SpinUDHalfOp(SpinUDHalfOp const & rhs) : SpinHalfOp() {
    this->index = rhs.index;
    this->plusF = rhs.plusF;
  }
  virtual ~SpinUDHalfOp() {}
  virtual std::string toString() const;
  virtual void herm() { plusF ^= 1; };
  virtual SpinHalfState operator*(SpinHalfBaseState const & rhs) const;
  virtual SpinHalfState operator*(SpinHalfState const & rhs) const;
};

//-------------------------------------------------------------------SpinMonomial--------

class SpinHalfMonomial : public MixMonomial<SpinHalfOp> {
 public:
  SpinHalfMonomial() : MixMonomial<SpinHalfOp>() {}
  SpinHalfMonomial(SpinHalfOp & Op) : MixMonomial<SpinHalfOp>(Op) {}
  SpinHalfMonomial(SpinHalfOp * OpPtr) : MixMonomial(OpPtr) {}
  //SpinHalfMonomial(Monomial<SpinHalfOp> const & rhs) : Monomial<SpinHalfOp>(rhs) {}
  ~SpinHalfMonomial() {}
  SpinHalfState operator*(SpinHalfBaseState const & rhs) const;
  SpinHalfState operator*(SpinHalfState const & rhs) const;
};

//------------------------------------------------------------------SpinPolynomial-------

class SpinHalfPolynomial : public Polynomial<SpinHalfOp> {};

#endif  //ORI_SDP_GS_SPINOPERATORS1D_HPP
