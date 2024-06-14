/*
  Jiazheng Sun
  Updated: Jun 13, 2024

  Class Implementations:
  DoubleDenseMatrix
*/

#ifndef LA_DENSE_NONTEM_CPP
#define LA_DENSE_NONTEM_CPP

#include <cblas.h>
#include <lapacke.h>

#include "./dense.hpp"

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
              ncols,
              rhs.data.data(),
              rhs.ncols,
              0.0,
              ans.data.data(),
              ans.ncols);
  return ans;
}

#endif  //LA_DENSE_NONTEM_CPP
