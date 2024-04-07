/*
  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_XXZ_HPP
#define ORI_SDP_GS_HAMILTONIANS_XXZ_HPP

#include "../Basics/hamiltonians.hpp"
#include "../Spin_1D/spinOperators1D.hpp"
#include "../Spin_1D/spinStates1D.hpp"

//------------------------------------------------------------XXZSparseHamiltonian-------

class XXZSparseHamiltonian
    : public SparseHamiltonian<SpinHalfPolynomial, SpinHalfBaseState> {
 private:
  size_t sites;
  double Jz;

 public:
  XXZSparseHamiltonian() : SparseHamiltonian<SpinHalfPolynomial, SpinHalfBaseState>() {}
  XXZSparseHamiltonian(SpinHalfPolynomial poly, size_t dim) :
      SparseHamiltonian<SpinHalfPolynomial, SpinHalfBaseState>(poly, dim) {}
  XXZSparseHamiltonian(SpinHalfPolynomial poly, size_t sites, double Jz) :
      SparseHamiltonian(poly, std::pow(2, sites)), sites(sites), Jz(Jz) {}
  ~XXZSparseHamiltonian() {}
};

//-------------------------------------------------------------XXZFullHamiltonian--------

class XXZFullHamiltonian : public FullHamiltonian<SpinHalfPolynomial, SpinHalfBaseState> {
 private:
  size_t sites;
  double Jz;

 public:
  XXZFullHamiltonian() : FullHamiltonian<SpinHalfPolynomial, SpinHalfBaseState>() {}
  XXZFullHamiltonian(SpinHalfPolynomial poly, size_t dim) :
      FullHamiltonian<SpinHalfPolynomial, SpinHalfBaseState>(poly, dim) {}
  XXZFullHamiltonian(SpinHalfPolynomial poly, size_t sites, double Jz) :
      FullHamiltonian(poly, std::pow(2, sites)), sites(sites), Jz(Jz) {}
  ~XXZFullHamiltonian() {}
  void initPoly(SpinHalfPolynomial poly) { this->poly = poly; };
  size_t getSites() const { return sites; }
  size_t getJz() const { return Jz; }
};

//----------------------------------------------------------------Other Functions--------

SpinHalfPolynomial makePoly(size_t sites, double Jz);

#endif  //ORI_SDP_GS_HAMILTONIANS_XXZ_HPP
