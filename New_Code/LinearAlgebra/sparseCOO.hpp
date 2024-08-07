/*
  Jiazheng Sun
  Updated: Aug 7, 2024
  
  Define tools for Coordinate Format (COO) sparse matrix operations.
*/

#ifndef LA_SPARSE_MATRIX_COO_HPP
#define LA_SPARSE_MATRIX_COO_HPP

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

//----------------------------------------------------------COOMatrix<DataType>----------

template<typename DataType>
class COOMatrix {
 protected:
  size_t nrows;                //Number of rows
  size_t ncols;                //Number of columns
  size_t nnz;                  //Number of non-zero elements
  std::vector<size_t> rows;    //Row indices of non-zero elements
  std::vector<size_t> cols;    //Column indices of non-zero elements
  std::vector<DataType> data;  //Values of non-zero elements

 public:
  COOMatrix() : nrows(0), ncols(0), nnz(0), rows(), cols(), data() {}
  COOMatrix(const size_t nrows, const size_t ncols) :
      nrows(nrows), ncols(ncols), nnz(0), rows(), cols(), data() {}
  COOMatrix(const COOMatrix<DataType> & rhs);
  virtual ~COOMatrix() {}
  /*Get information of the matrix.*/
  size_t getNrows() const { return nrows; }
  size_t getNcols() const { return ncols; }
  size_t getNnz() const { return nnz; }
  std::vector<size_t> getRows() const { return rows; }
  std::vector<size_t> getCols() const { return cols; }
  std::vector<DataType> getAllData() const { return data; }
  std::string toString() const;
  virtual std::string element_toString(DataType element) const = 0;
  /*Modify the matrix.*/
  void addData(size_t rowId, size_t colId, DataType newData);
};

//-----------------------------------------------------------ComplexCOOMatrix------------

class ComplexCOOMatrix : public COOMatrix<std::complex<double> > {
 public:
  ComplexCOOMatrix() : COOMatrix<std::complex<double> >() {}
  ComplexCOOMatrix(const size_t nrows, const size_t ncols) :
      COOMatrix<std::complex<double> >(nrows, ncols) {}
  ComplexCOOMatrix(const COOMatrix<std::complex<double> > & rhs) :
      COOMatrix<std::complex<double> >(rhs) {}
  virtual ~ComplexCOOMatrix() {}
  /*Get information of the matrix.*/
  virtual std::string element_toString(std::complex<double> element) const;
};

#endif  //LA_SPARSE_MATRIX_COO_HPP
