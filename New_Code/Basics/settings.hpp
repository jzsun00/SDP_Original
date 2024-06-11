/*
  Jiazheng Sun
  Updated: Jun 10, 2024

  Class:
  COOMatrix, CSRMatrix
  
  Function:
  string complex_toString(const complex<double> & num);
  string complexVector_toString(const vector<complex<double> > & vec);
  string complexMatrix_toString(const vector<vector<complex<double> > > & matrix);

  Define general settings that can be used for the entire project.
*/

#ifndef ORI_SDP_GS_SETTINGS_HPP
#define ORI_SDP_GS_SETTINGS_HPP

#include <complex>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

using std::complex;
using std::vector;

//------------------------------------------------------------Constant Definitions-------

/*Double and complex values less than ERROR will be considered zero in computations.*/
const double ERROR = std::pow(10, -12);
/*Default output precision for complex numbers.*/
const size_t COMPLEX_PRECISION = 4;

//--------------------------------------------------------------Sparse Matrices----------

template<typename DataType>
class COOMatrix {
 public:
  COOMatrix() {}
  ~COOMatrix() {}
};

template<typename DataType>
class CSRMatrix {};

//----------------------------------------------------------------Output Tools-----------

/*Convert a single complex number to std::string.*/
std::string complex_toString(const complex<double> & num);

/*Convert an std::vector of complex numbers to std::string.*/
std::string complexVector_toString(const vector<complex<double> > & vec);

/*Convert a dense matrix of complex numbers to std::string.*/
std::string complexMatrix_toString(const vector<vector<complex<double> > > & matrix);

#endif  //ORI_SDP_GS_SETTINGS_HPP
