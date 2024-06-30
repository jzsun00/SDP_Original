/*
  Jiazheng Sun
  Updated: Jun 27, 2024

  Class Implementations:
  COOMatrix<DataType>
  CSRMatrix<DataType>
*/

#ifndef ORI_SDP_GS_SETTINGS_TEM_CPP
#define ORI_SDP_GS_SETTINGS_TEM_CPP

#include "./settings.hpp"

//-----------------------------------------------------------------COOMatrix-------------

template<typename DataType>
void COOMatrix<DataType>::addData(size_t rowId, size_t colId, DataType newData) {
  rows.push_back(rowId);
  cols.push_back(colId);
  data.push_back(newData);
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

#endif  //ORI_SDP_GS_SETTINGS_TEM_CPP
