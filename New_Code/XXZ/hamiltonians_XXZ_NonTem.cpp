/*
  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
#define ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP

#include "./hamiltonians_XXZ.hpp"

SpinHalfPolynomial XXZFullHamiltonian::makePoly(size_t sites, double Jz) {
  SpinHalfPolynomial ans;
  for (size_t i = 0; i < sites - 1; i++) {
    SpinZHalfOp * SzPtr = new SpinZHalfOp(i);
    OpRecords.push_back(SzPtr);
    SpinZHalfOp * SzPtrN = new SpinZHalfOp(i + 1);
    OpRecords.push_back(SzPtrN);
    SpinUDHalfOp * SuPtr = new SpinUDHalfOp(i, true);
    OpRecords.push_back(SuPtr);
    SpinUDHalfOp * SdPtr = new SpinUDHalfOp(i, false);
    OpRecords.push_back(SdPtr);
    SpinUDHalfOp * SuPtrN = new SpinUDHalfOp(i + 1, true);
    OpRecords.push_back(SuPtrN);
    SpinUDHalfOp * SdPtrN = new SpinUDHalfOp(i + 1, false);
    OpRecords.push_back(SdPtrN);
    SpinHalfMonomial MNud(SuPtr);
    MNud *= SdPtrN;
    SpinHalfMonomial MNdu(SdPtr);
    MNdu *= SuPtrN;
    SpinHalfMonomial MNz(SzPtr);
    MNz *= SzPtrN;
    ans += pair<complex<double>, SpinHalfMonomial>(complex<double>(0.5, 0), MNud);
    ans += pair<complex<double>, SpinHalfMonomial>(complex<double>(0.5, 0), MNdu);
    ans += pair<complex<double>, SpinHalfMonomial>(complex<double>(Jz, 0), MNz);
  }
  return ans;
}

#endif  //ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
