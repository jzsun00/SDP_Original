/*
  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_HUBBARD_HPP
#define ORI_SDP_GS_HAMILTONIANS_HUBBARD_HPP

#include "../Basics/hamiltonians.hpp"
#include "../HardCore/hardCoreOperators.hpp"
#include "../HardCore/hardCoreSubspaces.hpp"

//----------------------------------------------------------------Other Functions--------

HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > makePoly_Hubbard1D(size_t sites,
                                                                             double V);

#endif  //ORI_SDP_GS_HAMILTONIANS_HUBBARD_HPP
