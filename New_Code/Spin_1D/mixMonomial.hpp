/*
  Jiazheng Sun
  Updated: Mar 12, 2024

  Define monomial with mixed types of operators.
  Designed to fit spin systems, where we have Sz, S+ and S- operators.
*/

#ifndef ORI_SDP_GS_MIXMONOMIAL_HPP
#define ORI_SDP_GS_MIXMONOMIAL_HPP

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using std::complex;
using std::pair;
using std::unique_ptr;
using std::vector;

template<typename OpType>
class MixMonomial {
 protected:
  vector<OpType *> ExprPtr;
  //std::unique_ptr<int> u3(new int);
 public:
  MixMonomial() : ExprPtr() {}
  MixMonomial(OpType & Op) : ExprPtr() { ExprPtr.push_back(&Op); }
  MixMonomial(OpType * OpPtr) : ExprPtr(1, OpPtr) {}
  ~MixMonomial() {
    for (typename vector<OpType *>::iterator it = ExprPtr.begin(); it != ExprPtr.end();
         ++it) {
      delete (*it);
    }
  }
  size_t getSize() const { return ExprPtr.size(); }
  std::string toString() const;
  OpType operator[](size_t n) const { return *(ExprPtr[n]); }
  MixMonomial & operator*=(OpType & toAdd);
  void addOp(OpType * OpPtr) { ExprPtr.push_back(OpPtr); }
  void herm();
};

template<typename OpType>
std::string MixMonomial<OpType>::toString() const {
  std::string ans = "";
  for (typename vector<OpType *>::const_iterator it = ExprPtr.begin();
       it != ExprPtr.end();
       ++it) {
    ans += (*it)->toString();
  }
  return ans;
}

template<typename OpType>
MixMonomial<OpType> & MixMonomial<OpType>::operator*=(OpType & toAdd) {
  ExprPtr.push_back(&toAdd);
  return *this;
}

template<typename OpType>
void MixMonomial<OpType>::herm() {
  std::reverse(ExprPtr.begin(), ExprPtr.end());
  for (typename vector<OpType *>::iterator it = ExprPtr.begin(); it != ExprPtr.end();
       ++it) {
    (*it)->herm();
  }
}

#endif
