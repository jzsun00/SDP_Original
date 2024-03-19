/*
  Jiazheng Sun
  Updated: Mar 18, 2024

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

#define ERROR std::pow(10, -12)
#define PRECISION 2

//-------------------------------------------------------------complex_toString()-------

std::string complex_toString(std::complex<double> num);

#endif  //ORI_SDP_GS_SETTINGS_HPP
