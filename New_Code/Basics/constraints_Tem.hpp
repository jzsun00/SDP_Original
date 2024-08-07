/*
  Jiazheng Sun
  Updated: Aug 6, 2024
  
  Class Implementations:
  ConsBaseSet<MonomialType, IndexType>
  ConsSet<MonomialType, IndexType>
  
  Function Implementations:
  void printMatrix(ConsSet & constraints, OpBasis & basis)
*/

#ifndef QM_CONSTRAINS_TEM_HPP
#define QM_CONSTRAINS_TEM_HPP

#include "./constraints.hpp"

//---------------------------------------ConsBaseSet<MonomialType, IndexType>------------

template<typename MonomialType, typename IndexType>
ConsBaseSet<MonomialType, IndexType> & ConsBaseSet<MonomialType, IndexType>::operator=(
    const ConsBaseSet<MonomialType, IndexType> & rhs) {
  if (this != &rhs) {
    this->start = rhs.start;
    this->end = rhs.end;
    this->order = rhs.order;
    this->BaseOpSet = rhs.BaseOpSet;
  }
  return *this;
}

//------------------------------------------ConsSet<MonomialType, IndexType>-------------

template<typename MonomialType, typename IndexType>
ConsSet<MonomialType, IndexType> & ConsSet<MonomialType, IndexType>::operator=(
    const ConsSet<MonomialType, IndexType> & rhs) {
  if (this != &rhs) {
    this->OpSet = rhs.OpSet;
  }
  return *this;
}

/*
template<typename MonomialType, typename IndexType>
Polynomial<MonomialType> ConsSet<MonomialType, IndexType>::getIJPoly(size_t i, size_t j) {
  MonomialType mnI = OpSet[i];
  mnI.herm();
  MonomialType mnJ = OpSet[j];
  mnI *= mnJ;
  Polynomial<MonomialType> ans(mnI);
  return ans;
}
*/

//-------------------------------------------------------Other Functions-----------------

template<typename MonomialType, typename IndexType>
void printMatrix(ConsSet<MonomialType, IndexType> & constrains,
                 OpBasis<MonomialType, IndexType> & basis) {
  size_t matrixNum = basis.getLength();
  size_t matrixSize = constrains.getLength();
  std::vector<std::vector<std::vector<std::complex<double> > > > matrices(
      matrixNum,
      std::vector<std::vector<std::complex<double> > >(
          matrixSize, std::vector<std::complex<double> >(matrixSize)));
  for (size_t i = 0; i < matrixSize; i++) {
    for (size_t j = 0; j < matrixSize; j++) {
      Polynomial<MonomialType> polyIJ = constrains.getIJPoly(i, j);
      std::vector<std::complex<double> > entryIJ = basis.projPoly(polyIJ);
      for (size_t k = 0; k < matrixNum; k++) {
        matrices[k][i][j] = entryIJ[k];
      }
    }
  }
}

#endif  //QM_CONSTRAINS_TEM_HPP
