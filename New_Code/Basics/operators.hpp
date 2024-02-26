/*
  Jiazheng Sun
  Updated: Feb 26, 2024
  Define ladder operators, monomials and polynomials.
  That serves as fundamental definitions of operators in the second quantization form,
  Fermi or Boson operators should inherit these classes and add related algebra.
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

class Monomial {
 protected:
  vector<LadderOp> Expr;

 public:
  /*Construct a monomial with one ladder operator or copy another,
   default constructor use an empty vector.*/
  Monomial() : Expr() {}
  Monomial(LadderOp & Op) : Expr(1, Op) {}
  Monomial(Monomial const & rhs) { Expr = rhs.Expr; }
  ~Monomial() {}
  /*Get information of the monomial.*/
  size_t getSize() const { return Expr.size(); }
  std::string toString();
  /*Overload operators.*/
  Monomial & operator=(Monomial const & rhs);
  bool operator==(Monomial const & rhs) const { return Expr == rhs.Expr; }
  LadderOp operator[](size_t n) const { return Expr[n]; }
  Monomial & operator*=(Monomial const & rhs);
  Monomial & operator*=(LadderOp const & toAdd);
  void herm();
};

//---------------------------------------------------------------Polynomial-------------

class Polynomial {
 protected:
  vector<pair<complex<double>, Monomial> > Terms;

 public:
  typedef pair<complex<double>, Monomial> TermType;
  /*Construct a polynomial with one monomial or copy another,
    default constructor use an empty vector.*/
  Polynomial() : Terms() {}
  Polynomial(Monomial const & mn) : Terms(1) {
    Terms[0].first = complex<double>(1, 0);
    Terms[0].second = mn;
  }
  Polynomial(Polynomial const & rhs) { Terms = rhs.Terms; }
  ~Polynomial() {}
  /*Get information of the polynomial.*/
  size_t getSize() const { return Terms.size(); }
  vector<pair<complex<double>, Monomial> >::iterator getBegin() { return Terms.begin(); }
  vector<pair<complex<double>, Monomial> >::iterator getEnd() { return Terms.end(); }
  std::string toString();
  /*Overload operators.*/
  TermType operator[](size_t n) const { return Terms[n]; }
  Polynomial & operator=(Polynomial const & rhs);
  bool operator==(Polynomial const & rhs) const { return Terms == rhs.Terms; }
  // If the prefactor of rhs is less than 10^(-9), ignore that operation
  Polynomial & operator+=(Monomial const & rhs);
  Polynomial & operator+=(TermType const & rhs);
  Polynomial & operator+=(Polynomial & rhs);
  Polynomial & operator-=(Monomial const & rhs);
  Polynomial & operator-=(TermType & rhs);
  Polynomial & operator-=(Polynomial & rhs);
  Polynomial & operator*=(Monomial const & rhs);
  Polynomial & operator*=(TermType const & rhs);
  Polynomial & operator*=(Polynomial & rhs);
  void herm();
  void eraseZeros();

 private:
  /*Find the same monomial for += operation.
    Return the corresponding iterator if same monomial is found,
    otherwise return Terms.end().*/
  vector<pair<complex<double>, Monomial> >::iterator findSameMonomial(
      Monomial const & mn);
};

#endif  //ORI_SDP_GS_OPERATORS_HPP
