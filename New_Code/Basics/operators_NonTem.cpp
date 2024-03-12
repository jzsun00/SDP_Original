/*
  Jiazheng Sun
  Updated: Mar 11, 2024
  Implementations of methods in class:
  LadderOp.
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

#endif  //ORI_SDP_GS_OPERATORS_NONTEM_CPP
