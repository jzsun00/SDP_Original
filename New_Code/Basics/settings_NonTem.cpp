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
  std::string ans;
  ans += std::to_string(num.real());
  ans += " + i";
  ans += std::to_string(num.imag());
  return ans;
}

#endif  //ORI_SDP_GS_SETTINGS_NONTEM_CPP
