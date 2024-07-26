/*
  Jiazheng Sun
  Updated: Jul 26, 2024

  Implementations of methods in class:
  OpBasis.
 */

#ifndef ORI_SDP_GS_SUBSPACES_TEM_HPP
#define ORI_SDP_GS_SUBSPACES_TEM_HPP

#include "./subspaces.hpp"

//-----------------------------------------------------------------OpSubBasis------------

template<typename MonomialType, typename IndexType>
std::vector<std::complex<double> > OpSubBasis<MonomialType, IndexType>::projPoly(
    Polynomial<MonomialType> poly) {
  std::vector<std::complex<double> > ans(Basis.size());
  for (size_t index = 0; index < Basis.size(); index++) {
    MonomialType basisMn = Basis[index];
    for (typename std::vector<
             std::pair<std::complex<double>, MonomialType> >::const_iterator it =
             poly.getBegin();
         it != poly.getEnd();
         ++it) {
      if (it->second == basisMn) {
        ans[index] = it->first;
      }
    }
  }
  return ans;
}

//--------------------------------------------------------------------OpBasis------------

template<typename MonomialType, typename IndexType>
std::vector<std::complex<double> > OpBasis<MonomialType, IndexType>::projPoly(
    Polynomial<MonomialType> poly) {
  std::vector<std::complex<double> > ans(Basis.size());
  for (size_t index = 0; index < Basis.size(); index++) {
    MonomialType basisMn = Basis[index];
    for (typename std::vector<
             std::pair<std::complex<double>, MonomialType> >::const_iterator it =
             poly.getBegin();
         it != poly.getEnd();
         ++it) {
      if (it->second == basisMn) {
        ans[index] = it->first;
      }
    }
  }
  return ans;
}

#endif  //ORI_SDP_GS_SUBSPACES_TEM_HPP
