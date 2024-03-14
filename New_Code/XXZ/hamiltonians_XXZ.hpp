/*
  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_XXZ_HPP
#define ORI_SDP_GS_HAMILTONIANS_XXZ_HPP

#include "../Basics/hamiltonians.hpp"
#include "../Spin_1D/spinOperators1D.hpp"
#include "../Spin_1D/spinStates1D.hpp"

class XXZFullHamiltonian : public FullHamiltonian<SpinHalfPolynomial, SpinHalfBasis> {
 private:
  size_t sites;
  double Jz;
  vector<SpinHalfOp *> OpRecords;

 public:
  XXZFullHamiltonian() : FullHamiltonian<SpinHalfPolynomial, SpinHalfBasis>() {}
  XXZFullHamiltonian(SpinHalfPolynomial poly, size_t dim) :
      FullHamiltonian<SpinHalfPolynomial, SpinHalfBasis>(poly, dim) {}
  XXZFullHamiltonian(SpinHalfPolynomial poly, size_t sites, double Jz) :
      FullHamiltonian(poly, std::pow(2, sites)), sites(sites), Jz(Jz), OpRecords() {}
  ~XXZFullHamiltonian() {}
  SpinHalfPolynomial makePoly(size_t sites, double Jz);
};

#endif  //ORI_SDP_GS_HAMILTONIANS_XXZ_HPP
