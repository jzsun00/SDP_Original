/*
  Jiazheng Sun
  Updated: Jun 12, 2024

  Class:
  DenseMatrix<DataType>

  Define general settings that can be used for the entire project.
*/

#ifndef LA_DENSE_HPP
#define LA_DENSE_HPP

#include <complex>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using std::complex;
using std::pair;
using std::vector;

/*Regular dense matrix.*/
template<typename DataType>
class DenseMatrix {
 protected:
  size_t nrows;                    //Number of rows
  size_t ncols;                    //Number of columns
  vector<vector<DataType> > data;  //Values of elements
 public:
  DenseMatrix() : nrows(0), ncols(0), data() {}
  DenseMatrix(size_t nrows, size_t ncols);
  DenseMatrix(const DenseMatrix<DataType> & rhs) :
      nrows(rhs.nrows), ncols(rhs.ncols), data(rhs.data) {}
  ~DenseMatrix() {}
  /*Get information of the matrix.*/
  pair<size_t, size_t> getDim() const { return pair<size_t, size_t>(nrows, ncols); }
  vector<vector<DataType> > getMatrix() const { return data; }
  DataType getData(size_t i, size_t j) const { return data[i][j]; }
  virtual std::string toString() const = 0;
};

#endif  //LA_DENSE_HPP
