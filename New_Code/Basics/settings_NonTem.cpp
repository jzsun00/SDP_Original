/*
  Jiazheng Sun
  Updated: Mar 18, 2024

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

#endif  //ORI_SDP_GS_SETTINGS_NONTEM_CPP
