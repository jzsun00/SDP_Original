/*
  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
#define ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP

#include "./hamiltonians_XXZ.hpp"

SpinHalfPolynomial makePoly(size_t sites, double Jz) {
  SpinHalfPolynomial ans;
  for (size_t i = 0; i < sites - 1; i++) {
    SpinHalfOp Sz(i);
    SpinHalfOp SzN(i + 1);
    SpinHalfOp Su(i, true);
    SpinHalfOp Sd(i, false);
    SpinHalfOp SuN(i + 1, true);
    SpinHalfOp SdN(i + 1, false);
    SpinHalfMonomial MNud(Su);
    MNud *= SdN;
    SpinHalfMonomial MNdu(Sd);
    MNdu *= SuN;
    SpinHalfMonomial MNz(Sz);
    MNz *= SzN;
    ans += pair<complex<double>, SpinHalfMonomial>(complex<double>(0.5, 0), MNud);
    ans += pair<complex<double>, SpinHalfMonomial>(complex<double>(0.5, 0), MNdu);
    ans += pair<complex<double>, SpinHalfMonomial>(complex<double>(Jz, 0), MNz);
  }
  return ans;
}

#endif  //ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
