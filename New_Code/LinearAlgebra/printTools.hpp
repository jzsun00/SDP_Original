/*
  Jiazheng Sun
  Updated: Aug 8, 2024

  Class:
  
  Function:
  string complex_toString(const complex<double> & num);
  string intVector_toString(const vector<int> & vec);
  string doubleVector_toString(const vector<double> & vec);
  string complexVector_toString(const vector<complex<double> > & vec);
  string complexMatrix_toString(const vector<vector<complex<double> > > & matrix);

  Define print tools to convert vectors and matrices to string.
*/

#ifndef LA_PRINT_TOOLS_HPP
#define LA_PRINT_TOOLS_HPP

#include <complex>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

//---------------------------------------------------------Constant Definitions-------

namespace LA {

  /*Default output precision for complex numbers.*/
  constexpr size_t COMPLEX_PRECISION = 4;
  /*Numerical value for Pi.*/
  constexpr double PI = 3.1415926536;

}  // namespace LA

//-------------------------------------------------------------Output Tools-----------

namespace LA {

  /*Convert a single complex number to std::string.*/
  std::string complex_toString(const std::complex<double> & num);

  /*Convert an std::vector of integer numbers to std::string.*/
  std::string intVector_toString(const std::vector<int> & vec);

  /*Convert an std::vector of size_t numbers to std::string.*/
  std::string size_tVector_toString(const std::vector<size_t> & vec);

  /*Convert an std::vector of double numbers to std::string.*/
  std::string doubleVector_toString(const std::vector<double> & vec);

  /*Convert an std::vector of complex numbers to std::string.*/
  std::string complexVector_toString(const std::vector<std::complex<double> > & vec);

  /*Convert a dense matrix of complex numbers to std::string.*/
  std::string complexMatrix_toString(
      const std::vector<std::vector<std::complex<double> > > & matrix);

}  // namespace LA

#endif  //LA_PRINT_TOOLS_HPP
