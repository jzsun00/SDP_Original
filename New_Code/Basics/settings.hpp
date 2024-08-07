/*
  Jiazheng Sun
  Updated: Aug 6, 2024

  Class:
  
  Function:
  string complex_toString(const complex<double> & num);
  string intVector_toString(const vector<int> & vec);
  string doubleVector_toString(const vector<double> & vec);
  string complexVector_toString(const vector<complex<double> > & vec);
  string complexMatrix_toString(const vector<vector<complex<double> > > & matrix);

  Define general settings that can be used for the entire project.
*/

#ifndef QM_SETTINGS_HPP
#define QM_SETTINGS_HPP

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

//-----------------------------------------------------------Constant Definitions-------

/*Double values less than ERROR will be considered zero in computations.*/
constexpr double ERROR = 1e-12;
/*Default output precision for complex numbers.*/
constexpr size_t COMPLEX_PRECISION = 4;
/*Numerical value for Pi.*/
constexpr double PI = 3.1415926536;

//----------------------------------------------------------------Output Tools----------

/*Convert a single complex number to std::string.*/
std::string complex_toString(const std::complex<double> & num);

/*Convert an std::vector of integer numbers to std::string.*/
std::string intVector_toString(const std::vector<int> & vec);

/*Convert an std::vector of double numbers to std::string.*/
std::string doubleVector_toString(const std::vector<double> & vec);

/*Convert an std::vector of complex numbers to std::string.*/
std::string complexVector_toString(const std::vector<std::complex<double> > & vec);

/*Convert a dense matrix of complex numbers to std::string.*/
std::string complexMatrix_toString(
    const std::vector<std::vector<std::complex<double> > > & matrix);

#endif  //QM_SETTINGS_HPP
