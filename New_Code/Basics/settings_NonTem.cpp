/*
  Jiazheng Sun
  Updated: Aug 15, 2024
  
  Function Implementations:
  string complex_toString(const complex<double> & num);
  string intVector_toString(const vector<int> & vec);
  string sizetVector_toString(const vector<size_t> & vec);
  string doubleVector_toString(const vector<double> & vec);
  string complexVector_toString(const vector<complex<double> > & vec);
  string complexMatrix_toString(const vector<vector<complex<double> > > & matrix);
*/

#ifndef QM_SETTINGS_NONTEM_CPP
#define QM_SETTINGS_NONTEM_CPP

#include <sstream>

#include "./settings.hpp"

using std::complex;
using std::vector;

//----------------------------------------------------------------Output Tools-----------

/*Convert a single complex number to std::string.*/
std::string complex_toString(const complex<double> & num) {
  std::ostringstream oss;
  oss.precision(COMPLEX_PRECISION);
  oss << std::scientific;
  oss << num.real();
  if (num.imag() > 0) {
    oss << "+" << num.imag() << "i";
  }
  else if (std::abs(num.imag()) < ERROR) {
    oss << "+" << 0 << "i";
  }
  else {
    oss << "-" << -num.imag() << "i";
  }
  return oss.str();
}

/*Convert an std::vector of integer numbers to std::string.*/
std::string intVector_toString(const vector<int> & vec) {
  std::string ans = "[ ";
  size_t count = 1;
  const size_t len = vec.size();
  for (size_t i = 0; i < len; ++i) {
    ans += (std::to_string(count) + "\t");
    ans += std::to_string(vec[i]);
    if (i != len - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

/*Convert an std::vector of size_t numbers to std::string.*/
std::string sizetVector_toString(const std::vector<size_t> & vec) {
  std::string ans = "[ ";
  size_t count = 1;
  const size_t len = vec.size();
  for (size_t i = 0; i < len; ++i) {
    ans += (std::to_string(count) + "\t");
    ans += std::to_string(vec[i]);
    if (i != len - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

/*Convert an std::vector of double numbers to std::string.*/
std::string doubleVector_toString(const vector<double> & vec) {
  std::string ans = "[ ";
  size_t count = 1;
  const size_t len = vec.size();
  for (size_t i = 0; i < len; ++i) {
    ans += (std::to_string(count) + "\t");
    ans += std::to_string(vec[i]);
    if (i != len - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

/*Convert an std::vector of complex numbers to std::string.*/
std::string complexVector_toString(const vector<complex<double> > & vec) {
  std::string ans = "[ ";
  size_t count = 1;
  const size_t len = vec.size();
  for (size_t i = 0; i < len; ++i) {
    ans += (std::to_string(count) + "\t");
    ans += complex_toString(vec[i]);
    if (i != len - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

/*Convert a dense matrix of complex numbers to std::string.*/
std::string complexMatrix_toString(const vector<vector<complex<double> > > & matrix) {
  std::string ans = "[ ";
  size_t count = 1;
  const size_t lenRow = matrix.size();
  const size_t lenCol = matrix[0].size();
  for (size_t i = 0; i < lenRow; ++i) {
    ans += (std::to_string(count) + "\t");
    for (size_t j = 0; j < lenCol; ++j) {
      ans += complex_toString(matrix[i][j]);
      ans += "\t";
    }
    if (i != lenRow - 1) {
      ans += "\n";
    }
    else {
      ans += "  ]";
    }
    count++;
  }
  return ans;
}

#endif  //QM_SETTINGS_NONTEM_CPP
