/*
  Jiazheng Sun
  Updated: Mar 20, 2024
*/

#ifndef ORI_SDP_GS_SUBSPACES_HPP
#define ORI_SDP_GS_SUBSPACES_HPP

#include "./operators.hpp"

template<typename MonomialType>
class OpBasis {
 protected:
  vector<MonomialType> Basis;

 public:
  OpBasis() : Basis() {}
  OpBasis(OpBasis const & rhs) : Basis(rhs.Basis) {}
  ~OpBasis() {}
  size_t getSize() const { return Basis.size(); }
  MonomialType operator[](size_t n) const { return Basis[n]; }
};

#endif  //ORI_SDP_GS_SUBSPACES_HPP
