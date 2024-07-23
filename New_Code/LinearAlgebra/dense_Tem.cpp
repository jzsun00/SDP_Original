/*
  Jiazheng Sun
  Updated: Jul 22, 2024

  Class Implementations:
  DenseMatrix<DataType>
*/

#ifndef LA_DENSE_MATRIX_TEM_CPP
#define LA_DENSE_MATRIX_TEM_CPP

#include <cstddef>

#include "./dense.hpp"

//--------------------------------------------------------DenseMatrix<DataType>-------

template<typename DataType>
DataType DenseMatrix<DataType>::operator()(size_t rowIdx, size_t colIdx) const {
  return data[colIdx * nrows + rowIdx];
}

template<typename DataType>
std::string DenseMatrix<DataType>::toString() const {
  std::string ans = "[ ";
  size_t count = 1;
  for (size_t rowId = 0; rowId < nrows; ++rowId) {
    ans += (std::to_string(count) + "\t");
    for (size_t colId = 0; colId < ncols; ++colId) {
      ans += data_toString(data[colId * nrows + rowId]);
      ans += "\t";
    }
    if (rowId != nrows - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

template<typename DataType>
void DenseMatrix<DataType>::setData(size_t rowIdx, size_t colIdx, DataType value) {
  data[colIdx * nrows + rowIdx] = value;
}

template<typename DataType>
bool DenseMatrix<DataType>::operator==(const DenseMatrix<DataType> & rhs) const {
  if (nrows != rhs.nrows || ncols != rhs.ncols) {
    return false;
  }
  return data = rhs.data;
}

template<typename DataType>
DenseMatrix<DataType> & DenseMatrix<DataType>::operator=(
    const DenseMatrix<DataType> & rhs) {
  nrows = rhs.nrows;
  ncols = rhs.ncols;
  data = rhs.data;
  return *this;
}

template<typename DataType>
DenseMatrix<DataType> & DenseMatrix<DataType>::operator+=(
    const DenseMatrix<DataType> & rhs) {
  const size_t len = nrows * ncols;
  for (size_t i = 0; i < len; i++) {
    data[i] += rhs.data[i];
  }
  return *this;
}

template<typename DataType>
DenseMatrix<DataType> & DenseMatrix<DataType>::operator-=(
    const DenseMatrix<DataType> & rhs) {
  const size_t len = nrows * ncols;
  for (size_t i = 0; i < len; i++) {
    data[i] -= rhs.data[i];
  }
  return *this;
}

#endif  //LA_DENSE_MATRIX_TEM_CPP
