/*
  Jiazheng Sun
  Updated: Jun 10, 2024

  Class:
  Operator, Ladderop, SpinOp, Monomial, Polynomial.

  Define general operators, monomials and polynomials.
  Also define ladder operators and spin operators inherited from operators.
  That serves as fundamental definitions of operators in the second quantization form
  as well as for spin systems.
  Fermi and Boson operators should inherit LadderOp, Monomial, Polynomial and
  add related algebra.
  Spin operators should inherit SpinOp, Monomial, Polynomial and add algebra.
*/

#ifndef ORI_SDP_GS_OPERATORS_HPP
#define ORI_SDP_GS_OPERATORS_HPP

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "./settings.hpp"

using std::complex;
using std::pair;
using std::vector;

//---------------------------------------------------------------Operator---------------

template<typename IndexType>
class Operator {
 protected:
  IndexType index;

 public:
  /*Construct a general operator with specified index.
    Default constructor use random value index.*/
  Operator() {}
  Operator(IndexType index) : index(index) {}
  Operator(Operator const & rhs) : index(rhs.index) {}
  ~Operator() {}
  /*Get information of the operator.*/
  IndexType getIndex() const { return index; }
  virtual std::string toString() const = 0;
  virtual std::string indexToString() const = 0;
  /*Overload operators.*/
  virtual void herm() = 0;
};

//---------------------------------------------------------------LadderOp---------------

template<typename IndexType>
class LadderOp : public Operator<IndexType> {
 protected:
  bool isUnit;
  bool creatorF;

 public:
  /*Construct a ladder operator with specified index and creatorF.
    Default constructor use random value index and false.*/
  LadderOp() : Operator<IndexType>(), isUnit(false), creatorF(false) {}
  LadderOp(IndexType index, bool creatorF) :
      Operator<IndexType>(index), isUnit(false), creatorF(creatorF) {}
  LadderOp(bool isUnit) : Operator<IndexType>(), isUnit(isUnit) {}
  LadderOp(LadderOp const & rhs) :
      Operator<IndexType>(rhs), isUnit(rhs.isUnit), creatorF(rhs.creatorF) {}
  ~LadderOp() {}
  /*Get information of the ladder operator.*/
  bool getCreatorF() const { return creatorF; }
  bool getIsUnit() const { return isUnit; }
  virtual std::string toString() const;
  virtual std::string indexToString() const = 0;
  /*Overload operators.*/
  virtual void herm() { creatorF ^= 1; }
};

//------------------------------------------------------------------SpinOp--------------

template<typename IndexType>
class SpinOp : public Operator<IndexType> {
 protected:
  bool isZ;
  bool isPlus;

 public:
  /*Construct a spin operator with specified index,
    default constructor use INT_MIN and Sz.
    If passing in one int, it's Sz;
    if passing in one int and one bool, it's S+ or S-.*/
  SpinOp() : Operator<IndexType>(), isZ(true), isPlus(false) {}
  SpinOp(IndexType index) : Operator<IndexType>(index), isZ(true), isPlus(false) {}
  SpinOp(IndexType index, bool isPlus) :
      Operator<IndexType>(index), isZ(false), isPlus(isPlus) {}
  SpinOp(IndexType index, bool isZ, bool isPlus) :
      Operator<IndexType>(index), isZ(isZ), isPlus(isPlus) {}
  SpinOp(SpinOp const & rhs) :
      Operator<IndexType>(rhs), isZ(rhs.isZ), isPlus(rhs.isPlus) {}
  ~SpinOp() {}
  /*Get information of the ladder operator.*/
  bool getIsZ() const { return isZ; }
  bool getIsPlus() const { return isPlus; }
  virtual std::string toString() const;
  virtual std::string indexToString() const = 0;
  /*Overload operators.*/
  SpinOp & operator=(SpinOp const & rhs);
  bool operator==(SpinOp const & rhs) const;
  virtual void herm();
};

//---------------------------------------------------------------Monomial---------------

template<typename OpType>
class Monomial {
 protected:
  vector<OpType> Expr;

 public:
  /*Construct a monomial with one ladder operator or copy another.
    Default constructor use an empty vector.*/
  Monomial() : Expr() {}
  Monomial(OpType & Op) : Expr(1, Op) {}
  Monomial(vector<OpType> const & Expr) : Expr(Expr) {}
  Monomial(Monomial<OpType> const & rhs) { Expr = rhs.Expr; }
  ~Monomial() {}
  /*Get information of the monomial.*/
  size_t getSize() const { return Expr.size(); }
  std::string toString() const;
  /*Overload operators.*/
  Monomial<OpType> & operator=(Monomial<OpType> const & rhs);
  bool operator==(Monomial<OpType> const & rhs) const { return Expr == rhs.Expr; }
  bool operator!=(Monomial<OpType> const & rhs) const { return !(*this == rhs); }
  OpType operator[](size_t n) const { return Expr[n]; }
  Monomial<OpType> & operator*=(Monomial<OpType> const & rhs);
  Monomial<OpType> & operator*=(OpType const & rhs);
  void herm();
};

//---------------------------------------------------------------Polynomial-------------

template<typename MonomialType>
class Polynomial {
 protected:
  vector<pair<complex<double>, MonomialType> > Terms;

 public:
  typedef pair<complex<double>, MonomialType> TermType;
  /*Construct a polynomial with one monomial or copy another.
    Default constructor use an empty vector.*/
  Polynomial() : Terms() {}
  Polynomial(MonomialType const & mn) : Terms() {
    Terms.push_back(TermType(complex<double>(1, 0), mn));
  }
  Polynomial(complex<double> pref, MonomialType const & mn) : Terms() {
    Terms.push_back(TermType(pref, mn));
  }
  Polynomial(Polynomial const & rhs) { Terms = rhs.Terms; }
  ~Polynomial() {}
  /*Get information of the polynomial.*/
  size_t getSize() const { return Terms.size(); }
  typename vector<pair<complex<double>, MonomialType> >::const_iterator getBegin() const {
    return Terms.begin();
  }
  typename vector<pair<complex<double>, MonomialType> >::const_iterator getEnd() const {
    return Terms.end();
  }
  std::string toString() const;
  /*Overload operators.*/
  Polynomial & operator=(Polynomial const & rhs);
  bool operator==(Polynomial const & rhs) { return Terms == rhs.Terms; }
  // If the prefactor of rhs is less than 10^(-12), ignore += and -= operations
  Polynomial & operator+=(MonomialType const & rhs);
  Polynomial & operator+=(TermType const & rhs);
  Polynomial & operator+=(Polynomial const & rhs);
  Polynomial & operator-=(MonomialType const & rhs);
  Polynomial & operator-=(TermType const & rhs);
  Polynomial & operator-=(Polynomial const & rhs);
  Polynomial & operator*=(MonomialType const & rhs);
  Polynomial & operator*=(TermType const & rhs);
  Polynomial & operator*=(Polynomial const & rhs);
  Polynomial & operator*=(complex<double> rhs);
  void herm();
  void eraseZeros();

 protected:
  /*Find the same monomial for += operation.
    Return the corresponding iterator if same monomial is found,
    otherwise return Terms.end().*/
  typename vector<pair<complex<double>, MonomialType> >::iterator findSameMonomial(
      MonomialType const & mn);
};

#include "operators_Tem.cpp"

#endif  //ORI_SDP_GS_OPERATORS_HPP
