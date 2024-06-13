/*
  Jiazheng Sun
  Updated: Jun 11, 2024

  Implementations of classes:
  DoubleDenseMatrix

  Define linear algebra tools.
*/

#ifndef ORI_SDP_GS_LINEARALGEBRA_NONTEM_HPP
#define ORI_SDP_GS_LINEARALGEBRA_NONTEM_HPP

#include <string>

#include "./linearAlgebra.hpp"
extern "C" {
#include <cblas.h>
}

//----------------------------------------------------------DoubleDenseMatrix-------------

std::string DoubleDenseMatrix::toString() const {
  std::string ans = "[ ";
  size_t count = 1;
  for (size_t i = 0; i < nrows; ++i) {
    ans += (std::to_string(count) + "    ");
    for (size_t j = 0; j < ncols; ++j) {
      ans += std::to_string(data[i][j]);
      ans += "    ";
    }
    if (i != nrows - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

#endif  //ORI_SDP_GS_LINEARALGEBRA_NONTEM_HPP
