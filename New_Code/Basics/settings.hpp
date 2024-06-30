/*
  Jiazheng Sun
  Updated: Jun 27, 2024

  Class:
  COOMatrix<DataType>
  CSRMatrix<DataType>
  
  Function:
  string complex_toString(const complex<double> & num);
  string intVector_toString(const vector<int> & vec);
  string complexVector_toString(const vector<complex<double> > & vec);
  string complexMatrix_toString(const vector<vector<complex<double> > > & matrix);

  Define general settings that can be used for the entire project.
*/

#ifndef ORI_SDP_GS_SETTINGS_HPP
#define ORI_SDP_GS_SETTINGS_HPP

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <complex>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using std::complex;
using std::pair;
using std::vector;

//------------------------------------------------------------Constant Definitions-------

/*Double values less than ERROR will be considered zero in computations.*/
const double ERROR = std::pow(10, -12);
/*Default output precision for complex numbers.*/
const size_t COMPLEX_PRECISION = 4;
/*Numerical value for Pi.*/
const double PI = 3.1415926536;

//------------------------------------------------------------------Matrices-------------

/*Coordinate Format (COO) sparse matrix.*/
template<typename DataType>
class COOMatrix {
 private:
  size_t nrows;           //Number of rows
  size_t ncols;           //Number of columns
  size_t nnz;             //Number of non-zero elements
  vector<size_t> rows;    //Row indices of non-zero elements
  vector<size_t> cols;    //Column indices of non-zero elements
  vector<DataType> data;  //Values of non-zero elements

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
  pair<size_t, size_t> getDim() const { return pair<size_t, size_t>(nrows, ncols); }
  size_t getNnz() const { return nnz; }
  vector<size_t> getRows() const { return rows; }
  vector<size_t> getCols() const { return cols; }
  vector<DataType> getData() const { return data; }
  std::string toString() const = 0;
  /*Modify the matrix.*/
  void addData(size_t rowId, size_t colId, DataType newData);
};

/*Compressed Sparse Row (CSR) sparse matrix.*/
template<typename DataType>
class CSRMatrix {};

template<typename T>
int findIndex(const std::vector<T> & vec, const T & value);

//----------------------------------------------------------------Output Tools----------

/*Convert a single complex number to std::string.*/
std::string complex_toString(const complex<double> & num);

/*Convert an std::vector of integer numbers to std::string.*/
std::string intVector_toString(const vector<int> & vec);

/*Convert an std::vector of complex numbers to std::string.*/
std::string complexVector_toString(const vector<complex<double> > & vec);

/*Convert a dense matrix of complex numbers to std::string.*/
std::string complexMatrix_toString(const vector<vector<complex<double> > > & matrix);

#endif  //ORI_SDP_GS_SETTINGS_HPP
