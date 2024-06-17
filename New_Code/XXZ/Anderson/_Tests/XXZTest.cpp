/*
  Jiazheng Sun
  Updated: Jun 17, 2024

  Calculate Anderson bound of XXZ model ground state energy.
*/

#include <chrono>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include "../../hamiltonians_XXZ.hpp"
#include "../include/arcomp.h"
#include "../include/arlnsmat.h"
#include "../include/arlscomp.h"
#include "../matrices/complex/lcompsol.h"

using std::cout;
using std::endl;

int main() {
  /*Set parameters sites and Jz.*/
  size_t sites = 12;
  double Jz = 0;
  //int dim = std::pow(2, sites);
  //std::cout << "dim = " << dim << std::endl;

  /*Construct polynomial and basis.*/
  SpinHalfPolynomial1D poly = makePoly(sites, Jz);
  SpinHalfBasis1D basis(sites);
  basis.init(0);
  std::cout << "Basis construction complete!" << std::endl;
  //std::cout << "Basis:" << std::endl << basis.toString() << std::endl;
  size_t dim = basis.getSize();
  std::cout << "dim = " << dim << std::endl;

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

  // Defining what we need: the 5 lowest eigenvalues of A.
  ARluCompStdEig<double> dprob(5L, A, "SR");

  // Finding eigenvalues and eigenvectors.
  auto start_solve = std::chrono::high_resolution_clock::now();
  dprob.FindEigenvectors();
  auto end_solve = std::chrono::high_resolution_clock::now();
  auto duration_solve =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_solve - start_solve);
  cout << "\nSolving Running Time: " << duration_solve.count() << " ms" << endl;

  // Printing solution.
  Solution(A, dprob);

}  // main
