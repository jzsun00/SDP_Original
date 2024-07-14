/*
  Jiazheng Sun
  Updated: Jul 13, 2024

  Class Implementations:
  DoubleDenseMatrix
*/

#ifndef LA_DENSE_NONTEM_CPP
#define LA_DENSE_NONTEM_CPP

#include <cblas.h>
#include <lapacke.h>

#include <cstddef>
#include <cstdlib>

#include "./dense.hpp"
using std::vector;

//-----------------------------------------------------------DoubleDenseMatrix--------

std::string DoubleDenseMatrix::toString() const {
  std::string ans = "[ ";
  size_t count = 1;
  for (size_t rowId = 0; rowId < nrows; ++rowId) {
    ans += (std::to_string(count) + "    ");
    for (size_t colId = 0; colId < ncols; ++colId) {
      ans += std::to_string(data[colId * nrows + rowId]);
      ans += "    ";
    }
    if (rowId != nrows - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

void DoubleDenseMatrix::fillRandomNum() {
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = std::rand() % rand_max;
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

#endif  //LA_DENSE_NONTEM_CPP
