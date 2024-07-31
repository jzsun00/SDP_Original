/*
  Jiazheng Sun
  Updated: Jul 30, 2024
  
  Class:
  XXZSparseHamiltonian
  XXZSparseRealHamiltonian
  XXZFullHamiltonian
  
  Define Hamiltonian matrices for 1D XXZ model.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_XXZ_HPP
#define ORI_SDP_GS_HAMILTONIANS_XXZ_HPP

#include <cstddef>

#include "../Basics/hamiltonians_Tem.hpp"
#include "../Spin_1D/spinOperators1D.hpp"
#include "../Spin_1D/spinStates1D.hpp"
#include "omp.h"

//--------------------------------------------------------XXZSparseHamiltonian-----------

class XXZSparseHamiltonian
    : public SparseHamiltonian<SpinHalfPolynomial1D, SpinHalfBaseState1D> {
 private:
  size_t sites;
  double Jz;

 public:
  XXZSparseHamiltonian() :
      SparseHamiltonian<SpinHalfPolynomial1D, SpinHalfBaseState1D>() {}
  XXZSparseHamiltonian(SpinHalfPolynomial1D poly, size_t dim) :
      SparseHamiltonian<SpinHalfPolynomial1D, SpinHalfBaseState1D>(poly, dim) {}
  XXZSparseHamiltonian(SpinHalfPolynomial1D poly, size_t sites, double Jz) :
      SparseHamiltonian(poly, std::pow(2, sites)), sites(sites), Jz(Jz) {}
  ~XXZSparseHamiltonian() {}
  virtual void createMatrix(SpinHalfBasis1D & basis);
  void createSymMatrix(SpinHalfBasis1D & basis);
};

//------------------------------------------------------XXZSparseRealHamiltonian---------

class XXZSparseRealHamiltonian
    : public SparseRealHamiltonian<SpinHalfPolynomial1D, SpinHalfBaseState1D> {
 private:
  size_t sites;
  double Jz;

 public:
  XXZSparseRealHamiltonian() :
      SparseRealHamiltonian<SpinHalfPolynomial1D, SpinHalfBaseState1D>() {}
  XXZSparseRealHamiltonian(const SpinHalfPolynomial1D & poly, size_t dim) :
      SparseRealHamiltonian<SpinHalfPolynomial1D, SpinHalfBaseState1D>(poly, dim) {}
  XXZSparseRealHamiltonian(const SpinHalfPolynomial1D & poly, size_t sites, double Jz) :
      SparseRealHamiltonian(poly), sites(sites), Jz(Jz) {}
  virtual ~XXZSparseRealHamiltonian() {}
  /*Get information of the sparse Hamiltonian.*/
  int * getIrowData() { return irow.data(); }
  int * getPcolData() { return pcol.data(); }
  double * getNzValData() { return nzVal.data(); }
  /*Use the specified basis to create matrix.*/
  virtual void createMatrix(SpinHalfBasis1D & basis);
  void createFullBasisMatrix(const SpinHalfBasis1D & basis);
  void createRefSymMatrix(SpinHalfBasis1D & basis);
};

//--------------------------------------------------------XXZFullHamiltonian-------------

class XXZFullHamiltonian
    : public FullHamiltonian<SpinHalfPolynomial1D, SpinHalfBaseState1D> {
 private:
  size_t sites;
  double Jz;

 public:
  XXZFullHamiltonian() : FullHamiltonian<SpinHalfPolynomial1D, SpinHalfBaseState1D>() {}
  XXZFullHamiltonian(SpinHalfPolynomial1D poly, size_t dim) :
      FullHamiltonian<SpinHalfPolynomial1D, SpinHalfBaseState1D>(poly, dim) {}
  XXZFullHamiltonian(SpinHalfPolynomial1D poly, size_t sites, double Jz) :
      FullHamiltonian(poly, std::pow(2, sites)), sites(sites), Jz(Jz) {}
  ~XXZFullHamiltonian() {}
  void initPoly(SpinHalfPolynomial1D poly) { this->poly = poly; };
  size_t getSites() const { return sites; }
  size_t getJz() const { return Jz; }
};

//----------------------------------------------------------------Other Functions--------

SpinHalfPolynomial1D makePoly(size_t sites, double Jz);

SpinHalfState1D makeMidState(size_t sites, double Jz, const SpinHalfBaseState1D & rhs);

#endif  //ORI_SDP_GS_HAMILTONIANS_XXZ_HPP
