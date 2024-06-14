/*
  Jiazheng Sun
  Updated: Jun 13, 2024

  Class:
  DenseMatrix<DataType>
  DoubleDenseMatrix

  Define tools for dense matrix operations using BLAS and LAPACK.
*/

#ifndef LA_DENSE_HPP
#define LA_DENSE_HPP

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

using std::vector;

//---------------------------------------------------------------DenseMatrix----------

/*General dense matrix.
  Use column-major ordering for BLAS and LAPACK interface.*/
template<typename DataType>
class DenseMatrix {
 protected:
  size_t nrows;           //Number of rows
  size_t ncols;           //Number of columns
  vector<DataType> data;  //Values of elements
 public:
  DenseMatrix() : nrows(0), ncols(0), data() {}
  DenseMatrix(size_t nrows, size_t ncols) :
      nrows(nrows), ncols(ncols), data(nrows * ncols) {}
  DenseMatrix(const DenseMatrix<DataType> & rhs) :
      nrows(rhs.nrows), ncols(rhs.ncols), data(rhs.data) {}
  ~DenseMatrix() {}
  /*Get information of the matrix.*/
  size_t getNrows() const { return nrows; }
  size_t getNcols() const { return ncols; }
  vector<DataType> getAllData() const { return data; }
  DataType operator()(size_t rowIdx, size_t colIdx) const;
  virtual std::string toString() const = 0;
  /*Modify the matrix.*/
  void setData(size_t rowIdx, size_t colIdx, DataType value);
};

//-----------------------------------------------------------DoubleDenseMatrix--------

/*Double valued dense matrix.*/
class DoubleDenseMatrix : public DenseMatrix<double> {
 public:
  DoubleDenseMatrix() : DenseMatrix<double>() {}
  DoubleDenseMatrix(size_t nrows, size_t ncols) : DenseMatrix<double>(nrows, ncols) {}
  DoubleDenseMatrix(const DoubleDenseMatrix & rhs) : DenseMatrix<double>(rhs) {}
  ~DoubleDenseMatrix() {}
  /*Get information of the matrix.*/
  virtual std::string toString() const;
  DoubleDenseMatrix operator*(const DoubleDenseMatrix & rhs) const;
};

#endif  //LA_DENSE_HPP
