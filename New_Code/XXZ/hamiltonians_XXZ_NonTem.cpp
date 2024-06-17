/*
  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
#define ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP

#include "./hamiltonians_XXZ.hpp"

SpinHalfPolynomial1D makePoly(size_t sites, double Jz) {
  SpinHalfPolynomial1D ans;
  for (size_t i = 0; i < sites - 1; i++) {
    SpinHalfOp1D Sz(i);
    SpinHalfOp1D SzN(i + 1);
    SpinHalfOp1D Su(i, true);
    SpinHalfOp1D Sd(i, false);
    SpinHalfOp1D SuN(i + 1, true);
    SpinHalfOp1D SdN(i + 1, false);
    SpinHalfMonomial1D MNud(Su);
    MNud *= SdN;
    SpinHalfMonomial1D MNdu(Sd);
    MNdu *= SuN;
    SpinHalfMonomial1D MNz(Sz);
    MNz *= SzN;
    if (i == 0 || i == sites - 2) {
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.25, 0), MNud);
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.25, 0), MNdu);
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(Jz / 2, 0), MNz);
    }
    else {
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.5, 0), MNud);
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.5, 0), MNdu);
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(Jz, 0), MNz);
    }
  }
  return ans;
}

#endif  //ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
