/*
  Jiazheng Sun
  Updated: Jun 16, 2024

  Calculate Anderson bound of XXZ model ground state energy.
*/

#include <cstdio>
#include <string>
#include <vector>

#include "../../../XXZ/hamiltonians_XXZ.hpp"
#include "../../hamiltonians_XXZ.hpp"
#include "../include/arcomp.h"
#include "../include/arlnsmat.h"
#include "../include/arlscomp.h"
#include "../matrices/complex/lcompsol.h"

using std::complex;
using std::vector;

int main() {
  /*Set parameters sites and Jz.*/
  size_t sites = 7;
  double Jz = 0;
  int dim = std::pow(2, sites);
  std::cout << "dim = " << dim << std::endl;

  /*Construct polynomial and basis.*/
  SpinHalfPolynomial poly = makePoly(sites, Jz);
  SpinHalfBasis basis(sites);
  basis.init();
  std::cout << "Basis construction complete!" << std::endl;

  /*Consruct sparse Hamiltonian.*/
  XXZSparseHamiltonian ham(poly, sites, Jz);
  ham.createMatrix(basis);
  std::cout << "Hamiltonian construction complete!" << std::endl;
  //std::cout << "Full Basis:\n" << basis.toString() << std::endl;

  int nnz = ham.getNumNonZero();
  std::cout << "nnz = " << nnz << std::endl;
  int * irow = new int[nnz];
  //irow = ham.getIrow().data();
  //std::cout << "size(irow) = " << ham.getIrow().size() << std::endl;
  int * pcol = new int[dim + 1];
  //pcol = ham.getPcol().data();
  //std::cout << "size(pcol) = " << ham.getPcol().size() << std::endl;
  complex<double> * valA = new complex<double>[nnz];
  //valA = ham.getNzVal().data();
  //std::cout << "size(valA) = " << ham.getNzVal().size() << std::endl;

  vector<int> irowVec = ham.getIrow();
  vector<int> pcolVec = ham.getPcol();
  vector<complex<double> > valAVec = ham.getNzVal();

  for (int i = 0; i < nnz; i++) {
    irow[i] = irowVec[i];
    valA[i] = valAVec[i];
  }
  for (int i = 0; i < dim + 1; i++) {
    pcol[i] = pcolVec[i];
  }

  ARluNonSymMatrix<complex<double>, double> A(dim, nnz, valA, irow, pcol);

  // Defining what we need: the four eigenvectors of A with largest magnitude.

  ARluCompStdEig<double> dprob(5L, A, "SR");

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

}  // main
