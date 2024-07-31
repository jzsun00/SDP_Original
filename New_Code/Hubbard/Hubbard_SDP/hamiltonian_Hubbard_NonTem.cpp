/*
  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_HUBBARD_NONTEM_CPP
#define ORI_SDP_GS_HAMILTONIANS_HUBBARD_NONTEM_CPP

#include "./hamiltonian_Hubbard.hpp"

HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > makePoly_Hubbard1D(size_t sites,
                                                                             double Jz) {
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > ans;
  for (int i = 0; i < (int)sites - 1; i++) {
    //SpinHalfOp Sz(i);
    //SpinHalfOp SzN(i + 1);
    HardCore1DLadderOp Su(i, true);
    HardCore1DLadderOp Sd(i, false);
    HardCore1DLadderOp SuN(i + 1, true);
    HardCore1DLadderOp SdN(i + 1, false);
    HardCoreMonomial<HardCore1DLadderOp> MNud(Su);
    MNud *= SdN;
    HardCoreMonomial<HardCore1DLadderOp> MNdu(Sd);
    MNdu *= SuN;
    HardCoreMonomial<HardCore1DLadderOp> MNz(Su);
    MNz *= Sd;
    HardCoreMonomial<HardCore1DLadderOp> MNzN(SuN);
    MNzN *= SdN;
    HardCoreMonomial<HardCore1DLadderOp> MNzz(MNz);
    MNzz *= MNzN;
    HardCore1DLadderOp unitOp(true);
    HardCoreMonomial<HardCore1DLadderOp> unitMn(unitOp);

    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(0.5, 0), MNud);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(-0.5, 0), MNdu);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(Jz, 0), MNzz);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(-0.5 * Jz, 0), MNz);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(-0.5 * Jz, 0), MNzN);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(0.25 * Jz, 0), unitMn);
  }
  return ans;
}

#endif  //ORI_SDP_GS_HAMILTONIANS_HUBBARD_NONTEM_CPP
