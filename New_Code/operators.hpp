/*
  Jan 31, 2024  Jiazheng Sun

  Define ladder operators
*/

#ifndef ORI_SDP_GS_OPERATORS_HPP
#define ORI_SDP_GS_OPERATORS_HPP

#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <string>

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
  std::string toString() const {
    std::string ans = "a_";
    ans += std::to_string(index);
    if (creatorF) {
      ans += "+";
    }
    return ans;
  }
  /*Overload operators.*/
  LadderOp & operator=(LadderOp const & rhs) {
    index = rhs.index;
    creatorF = rhs.creatorF;
    return *this;
  }
  bool operator==(LadderOp const & rhs) const {
    return creatorF == rhs.creatorF && index == rhs.index;
  }
  bool operator<(LadderOp const & rhs) const {
    if (creatorF != rhs.creatorF) {
      throw std::invalid_argument("Two operators are not the same kind!\n");
    }
    else {
      return index < rhs.index;
    }
  }
  bool operator>(LadderOp const & rhs) const { return !(*this < rhs || *this == rhs); }
  void herm() { creatorF ^= 1; }
};

class Monomial {
 protected:
  std::vector<LadderOp> Expr;
  
 public:
  /*Construct a monomial with one ladder operator or copy another,
   default constructor use an empty vector.*/
  Monomial(): Expr() {}
  Monomial(LadderOp & Op): Expr(1, Op) {}
  Monomial(Monomial const & rhs) {
    Expr = rhs.Expr;
  }
  ~Monomial() {}
  /*Get information of the monomial.*/
  size_t getSize() const {
    return Expr.size();
  }
  std::string toString() {
    std::string ans = "";
    for (std::vector<LadderOp>::iterator it = Expr.begin(); it != Expr.end(); ++it) {
      ans += it->toString();
    }
    return ans;
  }
  /*Overload operators.*/
  Monomial & operator=(Monomial const & rhs) {
    Expr = rhs.Expr;
    return *this;
  }
  bool operator==(Monomial const & rhs) const {
    return Expr == rhs.Expr;
  }
  LadderOp operator[](size_t n) const {
    return Expr[n];
  }
  Monomial & operator*=(Monomial const & rhs) {
    Expr.insert(Expr.end(), rhs.Expr.begin(), rhs.Expr.end());
    return *this;
  }
  Monomial & operator*=(LadderOp const & toAdd) {
    Expr.push_back(toAdd);
    return *this;
  }
  void herm() {
    std::reverse(Expr.begin(), Expr.end());
    for (std::vector<LadderOp>::iterator it = Expr.begin(); it != Expr.end(); ++it) {
      it->herm();
    }
  }
};

#endif
