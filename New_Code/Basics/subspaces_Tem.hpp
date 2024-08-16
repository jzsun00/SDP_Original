/*
  Jiazheng Sun
  Updated: Aug 15, 2024
  
  Class Implementations:
  OpSubBasis<MonomialType, IndexType>
  OpBasis<MonomialType, IndexType>
*/

#ifndef QM_SUBSPACES_TEM_HPP
#define QM_SUBSPACES_TEM_HPP

#include "./subspaces.hpp"

//-----------------------------------------OpSubBasis<MonomialType, IndexType>-----------

template<typename MonomialType, typename IndexType>
OpSubBasis<MonomialType, IndexType> & OpSubBasis<MonomialType, IndexType>::operator=(
    const OpSubBasis<MonomialType, IndexType> & rhs) {
  if (this != &rhs) {
    this->start = rhs.start;
    this->end = rhs.end;
    this->order = rhs.order;
    this->Basis = rhs.Basis;
  }
  return *this;
}

template<typename MonomialType, typename IndexType>
std::vector<std::complex<double> > OpSubBasis<MonomialType, IndexType>::projPoly(
    const Polynomial<MonomialType> & poly) const {
  const size_t len = Basis.size();
  std::vector<std::complex<double> > ans(len, std::complex<double>(0, 0));
  for (size_t index = 0; index < len; index++) {
    MonomialType basisMn = Basis[index];
    for (auto it = poly.getBegin(); it != poly.getEnd(); ++it) {
      if (it->second == basisMn) {
        ans[index] = it->first;
        break;
      }
    }
  }
  return ans;
}

//-------------------------------------------OpBasis<MonomialType, IndexType>------------

template<typename MonomialType, typename IndexType>
OpBasis<MonomialType, IndexType> & OpBasis<MonomialType, IndexType>::operator=(
    const OpBasis<MonomialType, IndexType> & rhs) {
  if (this != &rhs) {
    this->Basis = rhs.Basis;
  }
  return *this;
}

template<typename MonomialType, typename IndexType>
std::vector<std::complex<double> > OpBasis<MonomialType, IndexType>::projPoly(
    const Polynomial<MonomialType> & poly) const {
  const size_t len = Basis.size();
  std::vector<std::complex<double> > ans(len, std::complex<double>(0, 0));
  for (size_t index = 0; index < len; index++) {
    MonomialType basisMn = Basis[index];
    for (auto it = poly.getBegin(); it != poly.getEnd(); ++it) {
      if (it->second == basisMn) {
        ans[index] = it->first;
        break;
      }
    }
  }
  return ans;
}

#endif  //QM_SUBSPACES_TEM_HPP
