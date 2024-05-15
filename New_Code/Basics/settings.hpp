/*
  Jiazheng Sun
  Updated: May 14, 2024

  Class:
  
  Function:
  string complex_toString(complex<double> num);
  string complexVector_toString(vector<complex<double> > vec);
  string complexMatrix_toString(vector<vector<complex<double> > > matrix);

  Define general settings that can be used for the entire project.
*/

#ifndef ORI_SDP_GS_SETTINGS_HPP
#define ORI_SDP_GS_SETTINGS_HPP

#include <complex>
#include <sstream>
#include <string>
#include <vector>

//--------------------------------------------------------------Macro Definitions--------

/*Double and complex values below ERROR will be considered zero in computations.*/
#define ERROR std::pow(10, -12)
/*Default output precision for complex numbers.*/
#define COMPLEX_PRECISION 2

//----------------------------------------------------------------Output Tools-----------

/*Convert a single complex number to std::string.*/
std::string complex_toString(std::complex<double> num);

/*Convert an std::vector of complex numbers to std::string.*/
std::string complexVector_toString(std::vector<std::complex<double> > vec);

/*Convert a matrix of complex numbers to std::string.*/
std::string complexMatrix_toString(
    std::vector<std::vector<std::complex<double> > > matrix);

#endif  //ORI_SDP_GS_SETTINGS_HPP
