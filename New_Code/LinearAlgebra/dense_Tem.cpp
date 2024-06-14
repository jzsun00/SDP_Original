/*
  Jiazheng Sun
  Updated: Jun 13, 2024

  Class Implementations:
  DenseMatrix<DataType>
*/

#ifndef LA_DENSE_TEM_CPP
#define LA_DENSE_TEM_CPP

#include "./dense.hpp"

//---------------------------------------------------------------DenseMatrix----------

template<typename DataType>
DataType DenseMatrix<DataType>::operator()(size_t rowIdx, size_t colIdx) const {
  return data[colIdx * nrows + rowIdx];
}

template<typename DataType>
void DenseMatrix<DataType>::setData(size_t rowIdx, size_t colIdx, DataType value) {
  data[colIdx * nrows + rowIdx] = value;
}

#endif  //LA_DENSE_TEM_CPP
