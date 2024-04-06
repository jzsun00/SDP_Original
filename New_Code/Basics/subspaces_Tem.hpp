/*
  Jiazheng Sun
  Updated: Apr 6, 2024

  Implementations of methods in class:
  OpBasis.
 */

#ifndef ORI_SDP_GS_SUBSPACES_TEM_CPP
#define ORI_SDP_GS_SUBSPACES_TEM_CPP

#include "./subspaces.hpp"

//---------------------------------------------------------------OpBasis----------------

template<typename MonomialType, typename IndexType>
vector<complex<double> > OpSubBasis<MonomialType, IndexType>::projPoly(
    Polynomial<MonomialType> poly) {
  vector<complex<double> > ans(Basis.size());
  for (size_t index = 0; index < Basis.size(); index++) {
    MonomialType basisMn = Basis[index];
    for (typename vector<pair<complex<double>, MonomialType> >::const_iterator it =
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

#endif  //ORI_SDP_GS_SUBSPACES_TEM_CPP
