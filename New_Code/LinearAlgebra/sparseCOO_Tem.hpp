/*
  Jiazheng Sun
  Updated: Aug 8, 2024

  Class Implementations:
  COOMatrix<DataType>
*/

#ifndef LA_SPARSE_COO_TEM_HPP
#define LA_SPARSE_COO_TEM_HPP

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "./sparseCOO.hpp"

using std::complex;
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
std::string COOMatrix<DataType>::element_toString(DataType element) const {
  DataType copy = element;
  copy += element;
  return "";
}

template<typename DataType>
void COOMatrix<DataType>::addData(size_t rowId, size_t colId, DataType newData) {
  if (rowId >= nrows || colId >= ncols) {
    throw std::invalid_argument("ERROR: rowId or colId out of bounds!\n");
  }
  // Find the correct position to insert the new element
  auto it = std::lower_bound(rows.begin(), rows.end(), rowId);
  auto pos = std::distance(rows.begin(), it);
  while (pos < (long int)rows.size() && rows[pos] == rowId && cols[pos] < colId) {
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

template<typename DataType>
COOMatrix<DataType> & COOMatrix<DataType>::operator=(const COOMatrix<DataType> & rhs) {
  if (this != &rhs) {
    this->nrows = rhs.nrows;
    this->ncols = rhs.ncols;
    this->nnz = rhs.nnz;
    this->rows = rhs.rows;
    this->cols = rhs.cols;
    this->data = rhs.data;
  }
  return *this;
}

template<typename DataType>
COOMatrix<DataType> & COOMatrix<DataType>::operator+=(const COOMatrix<DataType> & rhs) {
  if (nrows != rhs.nrows || ncols != rhs.ncols) {
    throw std::invalid_argument("ERROR: Matrix dimensions must match for addition!\n");
  }
  COOMatrix<DataType> result(nrows, ncols);
  // Use a hash map to accumulate results
  std::unordered_map<size_t, DataType> resultMap;
  for (size_t i = 0; i < nnz; ++i) {
    size_t index = rows[i] * ncols + cols[i];
    resultMap[index] = data[i];
  }
  for (size_t i = 0; i < rhs.nnz; ++i) {
    size_t index = rhs.rows[i] * ncols + rhs.cols[i];
    resultMap[index] += rhs.data[i];
  }
  // Add accumulated results to the result matrix
  for (const auto & entry : resultMap) {
    size_t row = entry.first / ncols;
    size_t col = entry.first % ncols;
    DataType value = entry.second;
    result.addData(row, col, value);
  }
  *this = result;
  return *this;
}

template<typename DataType>
COOMatrix<DataType> & COOMatrix<DataType>::operator-=(const COOMatrix<DataType> & rhs) {
  if (nrows != rhs.nrows || ncols != rhs.ncols) {
    throw std::invalid_argument("ERROR: Matrix dimensions must match for subtraction!\n");
  }
  COOMatrix<DataType> result(nrows, ncols);
  // Use a hash map to accumulate results
  std::unordered_map<size_t, DataType> resultMap;
  for (size_t i = 0; i < nnz; ++i) {
    size_t index = rows[i] * ncols + cols[i];
    resultMap[index] = data[i];
  }
  for (size_t i = 0; i < rhs.nnz; ++i) {
    size_t index = rhs.rows[i] * ncols + rhs.cols[i];
    resultMap[index] -= rhs.data[i];
  }
  // Add accumulated results to the result matrix
  for (const auto & entry : resultMap) {
    size_t row = entry.first / ncols;
    size_t col = entry.first % ncols;
    DataType value = entry.second;
    result.addData(row, col, value);
  }
  *this = result;
  return *this;
}

template<typename DataType>
COOMatrix<DataType> & COOMatrix<DataType>::operator*=(const complex<double> pref) {
  for (size_t i = 0; i < nnz; i++) {
    data[i] *= pref;
  }
  return *this;
}

#endif  //LA_SPARSE_COO_TEM_HPP
