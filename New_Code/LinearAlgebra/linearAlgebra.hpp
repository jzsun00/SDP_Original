/*
  Jiazheng Sun
  Updated: Jun 11, 2024

  Class:
  DoubleDenseMatrix

  Define linear algebra tools.
*/

#ifndef ORI_SDP_GS_LINEARALGEBRA_HPP
#define ORI_SDP_GS_LINEARALGEBRA_HPP

#include <cstddef>

#include "settings.hpp"

//----------------------------------------------------------DoubleDenseMatrix-------------

/*Double valued dense matrix.*/
class DoubleDenseMatrix : public DenseMatrix<double> {
 public:
  DoubleDenseMatrix() : DenseMatrix<double>() {}
  DoubleDenseMatrix(size_t nrows, size_t ncols) : DenseMatrix<double>(nrows, ncols) {}
  DoubleDenseMatrix(const DoubleDenseMatrix & rhs) : DenseMatrix<double>(rhs) {}
  ~DoubleDenseMatrix() {}
  /*Get information of the matrix.*/
  virtual std::string toString() const;
  DoubleDenseMatrix operator*(DoubleDenseMatrix & rhs) const;
};

#endif  //ORI_SDP_GS_LINEARALGEBRA_HPP
