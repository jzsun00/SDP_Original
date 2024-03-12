/*
  Jiazheng Sun
  Updated: Mar 10, 2024

  Implementations of methods in class:
  FermiFockstate, FermiState, FermiBasis.
 */

#ifndef ORI_SDP_GS_SPINOPERATORS_CPP
#define ORI_SDP_GS_SPINOPERATORS_CPP

#include "spinOperators.hpp"

//------------------------------------------------------------------------SpinOp--------

template<typename SpinType>
SpinOp<SpinType> & SpinOp<SpinType>::operator=(SpinOp<SpinType> const & rhs) {
  index = rhs.index;
  type = rhs.type;
  return *this;
}

template<typename SpinType>
bool SpinOp<SpinType>::operator==(SpinOp<SpinType> const & rhs) const {
  return type == rhs.type && index == rhs.index;
}

template<typename SpinType>
bool SpinOp<SpinType>::operator<(SpinOp<SpinType> const & rhs) const {
  if (type != rhs.type) {
    throw std::invalid_argument("Two operators are not the same kind!\n");
  }
  else {
    return index < rhs.index;
  }
}

//----------------------------------------------------------------------SpinZOp---------

template<typename SpinType>
std::string SpinZOp<SpinType>::toString() const {
  std::string ans = "SZ_{";
  ans += std::to_string(this->index);
  ans += "}";
  return ans;
}

//----------------------------------------------------------------------SpinUDOp--------

template<typename SpinType>
std::string SpinUDOp<SpinType>::toString() const {
  std::string ans = "S";
  if (this->type == PLUS) {
    ans += "+";
  }
  if (this->type == MINUS) {
    ans += "-";
  }
  ans += std::to_string(this->index);
  ans += "}";
  return ans;
}

#endif
