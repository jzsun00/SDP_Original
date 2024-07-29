/*
  Jiazheng Sun
  Updated: Jul 28, 2024

  Calculate Anderson bound of XXZ model ground state energy.
*/

#include <chrono>
#include <cstddef>
#include <cstdio>

#include "../../hamiltonians_XXZ.hpp"
#include "../include/arcomp.h"
#include "../include/arlnsmat.h"
#include "../include/arlscomp.h"
#include "../matrices/complex/lcompsol.h"
#include "omp.h"

using std::complex;
using std::cout;
using std::endl;
using std::vector;

int main() {
  /*Set parameters sites and Jz.*/
  size_t sites = 4;
  double Jz = 0;
  cout << "Number of sites = " << sites << endl;
  cout << "Jz = " << Jz << endl << endl;

  omp_set_num_threads(12);

  /*Construct polynomial and basis.*/
  SpinHalfPolynomial1D poly = makePoly(sites, Jz);
  SpinHalfBasis1D * basis = new SpinHalfBasis1D(sites);
  auto start_basis_init = std::chrono::high_resolution_clock::now();
  basis->init(0);
  auto end_basis_init = std::chrono::high_resolution_clock::now();
  std::cout << "Basis construction complete!" << std::endl;
  //std::cout << "Basis:" << std::endl << basis.toString() << std::endl;
  size_t dim = basis->getSize();
  std::cout << "dim = " << dim << std::endl;

  /*Consruct sparse Hamiltonian.*/
  XXZSparseHamiltonian * ham = new XXZSparseHamiltonian(poly, sites, Jz);
  auto start_matrix_init = std::chrono::high_resolution_clock::now();
  ham->createMatrix(*basis);
  auto end_matrix_init = std::chrono::high_resolution_clock::now();
  std::cout << "Hamiltonian construction complete!" << std::endl;
  delete basis;
  //std::cout << "Full Basis:\n" << basis.toString() << std::endl;

  int nnz = ham->getNumNonZero();
  std::cout << "nnz = " << nnz << std::endl;
  //int * irow = new int[nnz];
  //irow = ham.getIrow().data();
  //std::cout << "size(irow) = " << ham.getIrow().size() << std::endl;
  //int * pcol = new int[dim + 1];
  //pcol = ham.getPcol().data();
  //std::cout << "size(pcol) = " << ham.getPcol().size() << std::endl;
  //complex<double> * valA = new complex<double>[nnz];
  //valA = ham.getNzVal().data();
  //std::cout << "size(valA) = " << ham.getNzVal().size() << std::endl;

  vector<int> irowVec = ham->getIrow();
  vector<int> pcolVec = ham->getPcol();
  vector<complex<double> > valAVec = ham->getNzVal();
  delete ham;
  //cout << "irow = " << intVector_toString(irowVec) << endl;
  //cout << "pcol = " << intVector_toString(pcolVec) << endl;
  //cout << "valA = " << complexVector_toString(valAVec) << endl;
  /*
  for (int i = 0; i < nnz; i++) {
    irow[i] = irowVec[i];
    valA[i] = valAVec[i];
  }
  for (int i = 0; i < dim + 1; i++) {
    pcol[i] = pcolVec[i];
  }
  */

  //ARluNonSymMatrix<complex<double>, double> A(dim, nnz, valA, irow, pcol);
  ARluNonSymMatrix<complex<double>, double> A(
      dim, nnz, valAVec.data(), irowVec.data(), pcolVec.data());

  // Defining what we need: the 3 lowest eigenvalues of A.
  ARluCompStdEig<double> dprob(5L, A, "SR");

  // Finding eigenvalues and eigenvectors.
  auto start_solve = std::chrono::high_resolution_clock::now();
  //omp_set_num_threads(4);
  dprob.FindEigenvectors();
  auto end_solve = std::chrono::high_resolution_clock::now();

  // Record time.
  auto duration_basis_init = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_basis_init - start_basis_init);
  auto duration_matrix_init = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_matrix_init - start_matrix_init);
  auto duration_solve =
      std::chrono::duration_cast<std::chrono::milliseconds>(end_solve - start_solve);
  cout << "\nInitiating Basis Running Time: " << duration_basis_init.count() << " ms"
       << endl;
  cout << "\nInitiating Matrix Running Time: " << duration_matrix_init.count() << " ms"
       << endl;
  cout << "\nSolving Running Time: " << duration_solve.count() << " ms" << endl;

  // Printing solution.
  Solution(A, dprob);
  double gs = dprob.Eigenvalue(0).real();
  for (size_t i = 1; i < 5; i++) {
    if (dprob.Eigenvalue(i).real() < gs) {
      gs = dprob.Eigenvalue(i).real();
    }
  }
  gs /= (sites - 2);
  cout << "\nGround State Energy = " << gs << endl;

}  // main
