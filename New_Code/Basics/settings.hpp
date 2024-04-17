/*
  Jiazheng Sun
  Updated: Apr 6, 2024

  Class:
  
  Function:
  std::string complex_toString(complex<double> num);

  Define general settings that can be used for the entire project.
*/

#ifndef ORI_SDP_GS_SETTINGS_HPP
#define ORI_SDP_GS_SETTINGS_HPP

#include <complex>
#include <sstream>
#include <string>
#include <vector>

#define ERROR std::pow(10, -12)
#define PRECISION 2

//-------------------------------------------------------------complex_toString()-------

std::string complex_toString(std::complex<double> num);

std::string complexVector_toString(std::vector<std::complex<double> > vec);

std::string complexMatrix_toString(
    std::vector<std::vector<std::complex<double> > > matrix);

#endif  //ORI_SDP_GS_SETTINGS_HPP
