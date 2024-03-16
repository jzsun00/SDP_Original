/*
  Jiazheng Sun
  Updated: Mar 16, 2024

  Implementations of methods in class:
  LadderOp, SpinOp.
 */

#ifndef ORI_SDP_GS_OPERATORS_NONTEM_CPP
#define ORI_SDP_GS_OPERATORS_NONTEM_CPP

#include "operators.hpp"

//---------------------------------------------------------------LadderOp---------------
std::string LadderOp::toString() const {
  std::string ans = "a_{";
  ans += std::to_string(index);
  ans += "}";
  if (creatorF) {
    ans += "{+}";
  }
  return ans;
}

LadderOp & LadderOp::operator=(LadderOp const & rhs) {
  index = rhs.index;
  creatorF = rhs.creatorF;
  return *this;
}

bool LadderOp::operator==(LadderOp const & rhs) const {
  return (creatorF == rhs.creatorF) && (index == rhs.index);
}

bool LadderOp::operator<(LadderOp const & rhs) const {
  if (creatorF != rhs.creatorF) {
    throw std::invalid_argument("Two operators are not the same kind!\n");
  }
  else {
    return index < rhs.index;
  }
}

//------------------------------------------------------------------SpinOp--------------

std::string SpinOp::toString() const {
  std::string ans = "S";
  if (isZ) {
    ans += "z";
  }
  else {
    if (isPlus) {
      ans += "+";
    }
    else {
      ans += "-";
    }
  }
  ans += "_{";
  ans += std::to_string(index);
  ans += "}";
  return ans;
}

SpinOp & SpinOp::operator=(SpinOp const & rhs) {
  index = rhs.index;
  isZ = rhs.isZ;
  isPlus = rhs.isPlus;
  return *this;
}

bool SpinOp::operator==(SpinOp const & rhs) const {
  return (index == rhs.index) && (isZ == rhs.isZ) && (isPlus == rhs.isPlus);
}

void SpinOp::herm() {
  if (isZ) {
    return;
  }
  else {
    isPlus ^= 1;
  }
}

//----------------------------------------------------------------Other Functions-------

std::string complex_toString(complex<double> num) {
  std::string ans;
  ans += std::to_string(num.real());
  ans += " + i";
  ans += std::to_string(num.imag());
  return ans;
}

#endif  //ORI_SDP_GS_OPERATORS_NONTEM_CPP
