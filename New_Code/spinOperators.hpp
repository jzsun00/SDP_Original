/*
  Jiazheng Sun
  Updated: Mar 10, 2024

  Define spin operators, monomials and polynomials for spin systems.
*/

#ifndef ORI_SDP_GS_SPINOPERATORS_HPP
#define ORI_SDP_GS_SPINOPERATORS_HPP

#include "./operators.hpp"
#include "./spinStates.hpp"

#define PLUS 0
#define MINUS 1
#define SZ 2

//------------------------------------------------------------------------SpinOp--------

template<typename SpinType>
class SpinOp {
 protected:
  int index;
  int type;

 public:
  /*Construct a spin operator with specified index and type,
   default constructor use INT_MIN and z-component.*/
  SpinOp() : index(std::numeric_limits<int>::min()), type(SZ) {}
  SpinOp(int index, int type) : index(index), type(type) {}
  SpinOp(SpinOp const & rhs) {
    index = rhs.index;
    type = rhs.type;
  }
  ~SpinOp() {}
  /*Get information of the spin operator.*/
  int getIndex() const { return index; }
  int getType() const { return type; }
  virtual std::string toString() const = 0;
  /*Overload operators.*/
  SpinOp & operator=(SpinOp const & rhs);
  bool operator==(SpinOp const & rhs) const;
  // Throw std::invalid_argument exception if operators are not the same kind
  bool operator<(SpinOp const & rhs) const;
  bool operator>(SpinOp const & rhs) const { return !(*this < rhs || *this == rhs); }
  virtual void herm() = 0;
  virtual SpinState<SpinType> operator*(SpinType & rhs) const;
  //virtual SpinState<SpinType> operator*(SpinState<SpinType> & rhs) const;
};

//----------------------------------------------------------------------SpinZOp---------

template<typename SpinType>
class SpinZOp : public SpinOp<SpinType> {
 public:
  SpinZOp() : SpinOp<SpinType>() {}
  SpinZOp(int index) : SpinOp<SpinType>(index, SZ) {}
  SpinZOp(SpinZOp const & rhs) : SpinOp<SpinType>(rhs) {}
  ~SpinZOp() {}
  virtual std::string toString() const;
  virtual void herm(){};
  virtual SpinState<SpinType> operator*(SpinType & rhs) const;
  //virtual SpinState<SpinType> operator*(SpinState<SpinType> & rhs) const;
};

//----------------------------------------------------------------------SpinUDOp--------

template<typename SpinType>
class SpinUDOp : public SpinOp<SpinType> {
 public:
  SpinUDOp() : SpinOp<SpinType>() {}
  SpinUDOp(int index, int type) : SpinOp<SpinType>(index, type) {}
  SpinUDOp(SpinUDOp const & rhs) : SpinOp<SpinType>(rhs) {}
  ~SpinUDOp() {}
  virtual std::string toString() const;
  virtual void herm();
  virtual SpinState<SpinType> operator*(SpinType & rhs) const;
  //virtual SpinState<SpinType> operator*(SpinState<SpinType> & rhs) const;
};

//-------------------------------------------------------------------SpinZHalfOp--------

class SpinZHalfOp : public SpinZOp<SpinHalfBaseState> {
 public:
  SpinZHalfOp() : SpinZOp() {}
  SpinZHalfOp(int index) : SpinZOp(index) {}
  SpinZHalfOp(SpinZOp const & rhs) : SpinZOp(rhs) {}
  ~SpinZHalfOp() {}
  virtual SpinState<SpinHalfBaseState> operator*(SpinHalfBaseState & rhs) const;
  //SpinHalfState operator*(SpinHalfState & rhs) const;
};

//-------------------------------------------------------------------SpinUDHalfOp-------

class SpinUDHalfOp : public SpinUDOp<SpinHalfBaseState> {
 public:
  SpinUDHalfOp() : SpinUDOp() {}
  SpinUDHalfOp(int index, int type) : SpinUDOp(index, type) {}
  SpinUDHalfOp(SpinUDOp const & rhs) : SpinUDOp(rhs) {}
  ~SpinUDHalfOp() {}
  virtual SpinState<SpinHalfBaseState> operator*(SpinHalfBaseState & rhs) const;
  //SpinHalfState operator*(SpinHalfState & rhs) const;
};

//-------------------------------------------------------------------SpinMonomial-------

template<typename SpinType>
class SpinMonomial : public Monomial<SpinOp<SpinType> > {};

//------------------------------------------------------------------SpinPolynomial------

template<typename SpinType>
class SpinPolynomial : public Polynomial<SpinOp<SpinType> > {};

#include "spinOperators.cpp"

#endif  //ORI_SDP_GS_SPINOPERATORS_HPP
