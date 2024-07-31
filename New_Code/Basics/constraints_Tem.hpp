/*
  Jiazheng Sun
  Updated: Jul 31, 2024
  
  Class Implementations:
  ConsBaseSet<MonomialType, IndexType>
  ConsSet<MonomialType, IndexType>
*/

#ifndef QM_CONSTRAINS_TEM_HPP
#define QM_CONSTRAINS_TEM_HPP

#include "./constraints.hpp"

//------------------------------------------ConsSet<MonomialType, IndexType>-------------
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
//-------------------------------------------------------------Other Functions-----------

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
