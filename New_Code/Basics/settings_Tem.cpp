/*
  Jiazheng Sun
  Updated: Jun 12, 2024

  Implementations of classes:
  DenseMatrix<DataType>
  COOMatrix<DataType>
  CSRMatrix<DataType>
*/

#ifndef ORI_SDP_GS_SETTINGS_TEM_CPP
#define ORI_SDP_GS_SETTINGS_TEM_CPP

#include "./settings.hpp"

//----------------------------------------------------------------DenseMatrix------------

template<typename DataType>
DenseMatrix<DataType>::DenseMatrix(size_t nrows, size_t ncols) :
    nrows(nrows),
    ncols(ncols),
    data(vector<vector<DataType> >(nrows, vector<DataType>(ncols))) {
}

//-----------------------------------------------------------------COOMatrix-------------

template<typename DataType>
void COOMatrix<DataType>::addData(size_t rowId, size_t colId, DataType newData) {
  rows.push_back(rowId);
  cols.push_back(colId);
  data.push_back(newData);
  nnz++;
}

#endif  //ORI_SDP_GS_SETTINGS_TEM_CPP
