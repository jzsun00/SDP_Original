/*
  Jiazheng Sun
  Updated: Jul 22, 2024

  Class:
  DenseMatrix<DataType>
  DoubleDenseMatrix

  Define tools for dense matrix operations using BLAS and LAPACK.
*/

#ifndef LA_DENSE_MATRIX_HPP
#define LA_DENSE_MATRIX_HPP

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

//--------------------------------------------------------------Constants-------------

//--------------------------------------------------------DenseMatrix<DataType>-------

/*General dense matrix.
  Use column-major ordering for BLAS and LAPACK interface.*/
template<typename DataType>
class DenseMatrix {
 protected:
  size_t nrows;                //Number of rows
  size_t ncols;                //Number of columns
  std::vector<DataType> data;  //Values of elements
 public:
  DenseMatrix() : nrows(0), ncols(0), data() {}
  DenseMatrix(const size_t nrows, const size_t ncols) :
      nrows(nrows), ncols(ncols), data(nrows * ncols) {}
  DenseMatrix(const DenseMatrix<DataType> & rhs) :
      nrows(rhs.nrows), ncols(rhs.ncols), data(rhs.data) {}
  virtual ~DenseMatrix() {}
  /*Get information of the matrix.*/
  size_t getNrows() const { return nrows; }
  size_t getNcols() const { return ncols; }
  std::vector<DataType> getAllData() const { return data; }
  DataType operator()(size_t rowIdx, size_t colIdx) const;
  std::string toString() const;
  virtual std::string data_toString(DataType element) const = 0;
  /*Modify the matrix.*/
  void setData(size_t rowIdx, size_t colIdx, DataType value);
  /*Operator overloading.*/
  bool operator==(const DenseMatrix<DataType> & rhs) const;
  DenseMatrix<DataType> & operator=(const DenseMatrix<DataType> & rhs);
  DenseMatrix<DataType> & operator+=(const DenseMatrix<DataType> & rhs);
  DenseMatrix<DataType> & operator-=(const DenseMatrix<DataType> & rhs);
};

//-----------------------------------------------------------DoubleDenseMatrix--------

/*Double valued dense matrix.*/
class DoubleDenseMatrix : public DenseMatrix<double> {
 public:
  DoubleDenseMatrix() : DenseMatrix<double>() {}
  DoubleDenseMatrix(const size_t nrows, const size_t ncols) :
      DenseMatrix<double>(nrows, ncols) {}
  DoubleDenseMatrix(const DoubleDenseMatrix & rhs) : DenseMatrix<double>(rhs) {}
  virtual ~DoubleDenseMatrix() {}
  /*Get information of the matrix.*/
  virtual std::string data_toString(double element) const;
  /*Modify the matrix.*/
  void fillRandomNum();
  void fillRandomNum(size_t maxNum);
  /*BLAS (OpenBLAS) tools.*/
  DoubleDenseMatrix operator*(const DoubleDenseMatrix & rhs) const;
  /*LAPACK tools.*/
  void solveEigen(std::vector<double> & real, std::vector<double> & imag);
};

#endif  //LA_DENSE_MATRIX_HPP
