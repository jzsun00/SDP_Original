/*
  Jiazheng Sun
  Updated: Jul 22, 2024

  Class Implementations:
  COOMatrix<DataType>
*/

#ifndef LA_SPARSE_COO_TEM_CPP
#define LA_SPARSE_COO_TEM_CPP

#include <stdexcept>

#include "./sparseCOO.hpp"

using std::vector;

//----------------------------------------------------------COOMatrix<DataType>----------

template<typename DataType>
COOMatrix<DataType>::COOMatrix(const COOMatrix<DataType> & rhs) :
    nrows(rhs.nrows),
    ncols(rhs.ncols),
    nnz(rhs.nnz),
    rows(rhs.rows),
    cols(rhs.cols),
    data(rhs.data) {
}

template<typename DataType>
void COOMatrix<DataType>::addData(size_t rowId, size_t colId, DataType newData) {
  if (rowId >= nrows || colId >= ncols) {
    throw std::invalid_argument("ERROR: rowId or colId out of bounds!\n");
  }
  // Find the correct position to insert the new element
  auto it = std::lower_bound(rows.begin(), rows.end(), rowId);
  auto pos = std::distance(rows.begin(), it);
  while (pos < rows.size() && rows[pos] == rowId && cols[pos] < colId) {
    ++pos;
  }
  rows.insert(rows.begin() + pos, rowId);
  cols.insert(cols.begin() + pos, colId);
  data.insert(data.begin() + pos, newData);
  nnz++;
}

template<typename T>
int findIndex(const std::vector<T> & vec, const T & value) {
  auto it = std::lower_bound(vec.begin(), vec.end(), value);

  // Check if the element is present in the vector
  if (it != vec.end() && *it == value) {
    return std::distance(vec.begin(), it);
  }
  else {
    return -1;  // Element not found
  }
}

#endif  //LA_SPARSE_COO_TEM_CPP
