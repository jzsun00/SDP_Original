/*
  Jiazheng Sun
  Updated: Mar 16, 2024

  Define ladder operators, monomials and polynomials.
  That serves as fundamental definitions of operators in the second quantization form,
  Fermi or Boson operators should inherit these classes and add related algebra.
  Spin operators can also inherit Monomial and Polynomial class.

  Implementations for LadderOp is in operators_Tem.cpp,
  for Monomial and Polynomial are in operators_NonTem.cpp.
  Only include _Tem.cpp at the end of this file.
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

#define ERROR std::pow(10, -12)

using std::complex;
using std::pair;
using std::vector;

//---------------------------------------------------------------LadderOp---------------

class LadderOp {
 protected:
  int index;
  bool creatorF;

 public:
  /*Construct a ladder operator with specified index and creatorF,
   default constructor use INT_MIN and false.*/
  LadderOp() : index(std::numeric_limits<int>::min()), creatorF(false) {}
  LadderOp(int i_idx, bool i_creatorF) : index(i_idx), creatorF(i_creatorF) {}
  LadderOp(LadderOp const & rhs) {
    index = rhs.index;
    creatorF = rhs.creatorF;
  }
  ~LadderOp() {}
  /*Get information of the ladder operator.*/
  int getIndex() const { return index; }
  bool getCreatorF() const { return creatorF; }
  std::string toString() const;
  /*Overload operators.*/
  LadderOp & operator=(LadderOp const & rhs);
  bool operator==(LadderOp const & rhs) const;
  // Throw std::invalid_argument exception if operators are not the same kind
  bool operator<(LadderOp const & rhs) const;
  bool operator>(LadderOp const & rhs) const { return !(*this < rhs || *this == rhs); }
  void herm() { creatorF ^= 1; }
};

//---------------------------------------------------------------Monomial---------------

template<typename OpType>
class Monomial {
 protected:
  vector<OpType> Expr;

 public:
  /*Construct a monomial with one ladder operator or copy another,
   default constructor use an empty vector.*/
  Monomial() : Expr() {}
  Monomial(OpType & Op) : Expr(1, Op) {}
  Monomial(Monomial const & rhs) { Expr = rhs.Expr; }
  ~Monomial() {}
  /*Get information of the monomial.*/
  size_t getSize() const { return Expr.size(); }
  std::string toString() const;
  /*Overload operators.*/
  Monomial & operator=(Monomial const & rhs);
  bool operator==(Monomial const & rhs) const { return Expr == rhs.Expr; }
  OpType operator[](size_t n) const { return Expr[n]; }
  Monomial & operator*=(Monomial const & rhs);
  Monomial & operator*=(OpType const & toAdd);
  void herm();
};

//---------------------------------------------------------------Polynomial-------------

template<typename MonomialType>
class Polynomial {
 protected:
  vector<pair<complex<double>, MonomialType> > Terms;

 public:
  typedef pair<complex<double>, MonomialType> TermType;
  /*Construct a polynomial with one monomial or copy another,
    default constructor use an empty vector.*/
  Polynomial() : Terms() {}
  Polynomial(MonomialType const & mn) : Terms() {
    Terms.push_back(TermType(complex<double>(1, 0), mn));
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
  std::string toString();
  /*Overload operators.*/
  Polynomial & operator=(Polynomial const & rhs);
  bool operator==(Polynomial const & rhs) const { return Terms == rhs.Terms; }
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
  void herm();
  void eraseZeros();

 protected:
  /*Find the same monomial for += operation.
    Return the corresponding iterator if same monomial is found,
    otherwise return Terms.end().*/
  virtual typename vector<pair<complex<double>, MonomialType> >::iterator
  findSameMonomial(MonomialType const & mn);
};

#include "operators_Tem.cpp"

#endif  //ORI_SDP_GS_OPERATORS_HPP
