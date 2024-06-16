/*
  Jiazheng Sun
  Updated: Jun 16, 2024

  Class:
  Operator<IndexType>
  LadderOp<IndexType>
  SpinOp<IndexType>
  Monomial<OpType>
  Polynomial<MonomialType>
  
  Define general operators, monomials and polynomials.
  Also define ladder operators and spin operators inherited from operators.
  That serves as fundamental definitions of operators in the second quantization
  form as well as for spin systems.
  Fermi and Boson operators should inherit LadderOp, Monomial, Polynomial and
  add related algebra.
  Spin operators should inherit SpinOp, Monomial, Polynomial and add algebra.
*/

#ifndef QM_OPERATORS_HPP
#define QM_OPERATORS_HPP

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <stdexcept>

#include "./settings.hpp"
#include "./settings_Tem.cpp"

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
  Operator(const Operator<IndexType> & rhs) : index(rhs.index) {}
  virtual ~Operator() {}
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
  virtual ~LadderOp() {}
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
  bool isUnit;

 public:
  /*Construct a spin operator with specified index,
    default constructor use random value index and Sz.
    If only passing in a bool value, it's identity operator;
    If only passing in an index, it's Sz with the specified index;
    if passing in one index and one bool, it's S+ or S-.*/
  SpinOp() : Operator<IndexType>(), isZ(true), isPlus(false), isUnit(false) {}
  SpinOp(bool isUnit) :
      Operator<IndexType>(), isZ(false), isPlus(false), isUnit(isUnit) {}
  SpinOp(IndexType index) :
      Operator<IndexType>(index), isZ(true), isPlus(false), isUnit(false) {}
  SpinOp(IndexType index, bool isPlus) :
      Operator<IndexType>(index), isZ(false), isPlus(isPlus), isUnit(false) {}
  SpinOp(IndexType index, bool isZ, bool isPlus) :
      Operator<IndexType>(index), isZ(isZ), isPlus(isPlus), isUnit(false) {}
  SpinOp(const SpinOp<IndexType> & rhs) :
      Operator<IndexType>(rhs), isZ(rhs.isZ), isPlus(rhs.isPlus), isUnit(rhs.isUnit) {}
  virtual ~SpinOp() {}
  /*Get information of the spin operator.*/
  bool getIsZ() const { return isZ; }
  bool getIsPlus() const { return isPlus; }
  bool getIsUnit() const { return isUnit; }
  virtual std::string toString() const;
  virtual std::string indexToString() const = 0;
  /*Overload operators.*/
  SpinOp<IndexType> & operator=(const SpinOp<IndexType> & rhs);
  bool operator==(const SpinOp & rhs) const;
  virtual void herm();
};

//---------------------------------------------------------------Monomial---------------

template<typename OpType>
class Monomial {
 protected:
  vector<OpType> Expr;

 public:
  /*Construct a monomial with one operator or a vector of operators.
    Default constructor use an empty vector.*/
  Monomial() : Expr() {}
  Monomial(OpType & Op) : Expr(1, Op) {}
  Monomial(const vector<OpType> & Expr) : Expr(Expr) {}
  Monomial(const Monomial<OpType> & rhs) : Expr(rhs.Expr) {}
  virtual ~Monomial() {}
  /*Get information of the monomial.*/
  size_t getSize() const { return Expr.size(); }
  std::string toString() const;
  /*Overload operators.*/
  Monomial<OpType> & operator=(const Monomial<OpType> & rhs);
  bool operator==(const Monomial<OpType> & rhs) const { return Expr == rhs.Expr; }
  bool operator!=(const Monomial<OpType> & rhs) const { return !(*this == rhs); }
  OpType operator[](size_t n) const { return Expr[n]; }
  Monomial<OpType> & operator*=(const Monomial<OpType> & rhs);
  Monomial<OpType> & operator*=(const OpType & rhs);
  void herm();
};

//---------------------------------------------------------------Polynomial-------------

template<typename MonomialType>
class Polynomial {
 protected:
  vector<pair<complex<double>, MonomialType> > Terms;

 public:
  typedef pair<complex<double>, MonomialType> TermType;
  /*Construct a polynomial with one monomial.
    Default constructor use an empty vector.*/
  Polynomial() : Terms() {}
  Polynomial(const MonomialType & mn) : Terms(1, TermType(complex<double>(1.0, 0), mn)) {}
  Polynomial(complex<double> pref, const MonomialType & mn) :
      Terms(1, TermType(pref, mn)) {}
  Polynomial(const Polynomial & rhs) : Terms(rhs.Terms) {}
  virtual ~Polynomial() {}
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

#endif  //QM_OPERATORS_HPP
