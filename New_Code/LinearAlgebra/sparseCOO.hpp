/*
  Jiazheng Sun
  Updated: Jul 13, 2024
  
  
*/

#ifndef LA_SPARSE_COO_HPP
#define LA_SPARSE_COO_HPP

#include <cstddef>
#include <string>
#include <vector>

//-----------------------------------------------------------------COOMatrix-------------

/*Coordinate Format (COO) sparse matrix.*/
template<typename DataType>
class COOMatrix {
 private:
  size_t nrows;                //Number of rows
  size_t ncols;                //Number of columns
  size_t nnz;                  //Number of non-zero elements
  std::vector<size_t> rows;    //Row indices of non-zero elements
  std::vector<size_t> cols;    //Column indices of non-zero elements
  std::vector<DataType> data;  //Values of non-zero elements

 public:
  COOMatrix() : nrows(0), ncols(0), nnz(0), rows(), cols(), data() {}
  COOMatrix(size_t nrows, size_t ncols) :
      nrows(nrows), ncols(ncols), nnz(0), rows(), cols(), data() {}
  COOMatrix(const COOMatrix<DataType> & rhs) :
      nrows(rhs.nrows),
      ncols(rhs.ncols),
      nnz(rhs.nnz),
      rows(rhs.rows),
      cols(rhs.cols),
      data(rhs.data) {}
  ~COOMatrix() {}
  /*Get information of the matrix.*/
  std::pair<size_t, size_t> getDim() const {
    return std::pair<size_t, size_t>(nrows, ncols);
  }
  size_t getNnz() const { return nnz; }
  std::vector<size_t> getRows() const { return rows; }
  std::vector<size_t> getCols() const { return cols; }
  std::vector<DataType> getData() const { return data; }
  std::string toString() const = 0;
  /*Modify the matrix.*/
  void addData(size_t rowId, size_t colId, DataType newData);
};

#endif  //LA_SPARSE_COO_HPP
