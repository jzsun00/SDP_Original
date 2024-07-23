/*
  Jiazheng Sun
  Updated: Jul 22, 2024

  Class Implementations:
  DoubleDenseMatrix
*/

#ifndef LA_DENSE_MATRIX_NONTEM_CPP
#define LA_DENSE_MATRIX_NONTEM_CPP

#include <cblas.h>
#include <lapacke.h>

#include "./dense.hpp"

//-----------------------------------------------------------DoubleDenseMatrix--------

std::string DoubleDenseMatrix::data_toString(double element) const {
  return std::to_string(element);
}

void DoubleDenseMatrix::fillRandomNum() {
  const size_t len = data.size();
  for (size_t i = 0; i < len; ++i) {
    data[i] = std::rand();
  }
}

void DoubleDenseMatrix::fillRandomNum(size_t maxNum) {
  const size_t len = data.size();
  for (size_t i = 0; i < len; ++i) {
    data[i] = std::rand() % maxNum;
  }
}

DoubleDenseMatrix DoubleDenseMatrix::operator*(const DoubleDenseMatrix & rhs) const {
  DoubleDenseMatrix ans(nrows, rhs.ncols);
  cblas_dgemm(CblasColMajor,
              CblasNoTrans,
              CblasNoTrans,
              nrows,
              rhs.ncols,
              ncols,
              1.0,
              data.data(),
              nrows,
              rhs.data.data(),
              rhs.nrows,
              0.0,
              ans.data.data(),
              ans.nrows);
  return ans;
}

void DoubleDenseMatrix::solveEigen(std::vector<double> & real,
                                   std::vector<double> & imag) {
  if (nrows != ncols) {
    throw std::invalid_argument("ERROR: Matrix must be square to solve eigenvalues!\n");
  }
  int n = nrows;
  int status = LAPACKE_dgeev(LAPACK_ROW_MAJOR,
                             'N',
                             'N',
                             n,
                             data.data(),
                             n,
                             real.data(),
                             imag.data(),
                             nullptr,
                             n,
                             nullptr,
                             n);
  if (status != 0) {
    throw std::invalid_argument("ERROR: LAPACKE_dgeev returns non-zero value " +
                                std::to_string(status) + " !\n");
  }
}

#endif  //LA_DENSE_MATRIX_NONTEM_CPP
