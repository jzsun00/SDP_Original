/*
  Jiazheng Sun
  Updated: Apr 6, 2024

  Implementations of functions:
  complex_toString().
*/

#ifndef ORI_SDP_GS_SETTINGS_NONTEM_CPP
#define ORI_SDP_GS_SETTINGS_NONTEM_CPP

#include "settings.hpp"

//-------------------------------------------------------------complex_toString()-------

std::string complex_toString(std::complex<double> num) {
  std::ostringstream oss;
  oss.precision(PRECISION);
  oss << std::fixed;
  oss << num.real();
  if (num.imag() >= 0) {
    oss << " + " << num.imag() << "i";
  }
  else {
    oss << " - " << -num.imag() << "i";
  }
  return oss.str();
}

std::string complexVector_toString(std::vector<std::complex<double> > vec) {
  std::string ans = "[ ";
  size_t count = 1;
  for (size_t i = 0; i < vec.size(); ++i) {
    ans += (std::to_string(count) + "    ");
    ans += complex_toString(vec[i]);
    if (i < vec.size() - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

std::string complexMatrix_toString(
    std::vector<std::vector<std::complex<double> > > matrix) {
  std::string ans = "[ ";
  size_t count = 1;
  for (size_t i = 0; i < matrix.size(); ++i) {
    ans += (std::to_string(count) + "    ");
    for (size_t j = 0; j < matrix.size(); ++j) {
      ans += complex_toString(matrix[i][j]);
      ans += "    ";
    }
    if (i < matrix.size() - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

#endif  //ORI_SDP_GS_SETTINGS_NONTEM_CPP
