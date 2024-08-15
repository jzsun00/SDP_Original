/*
  Jiazheng Sun
  Updated: Aug 15, 2024
  
  Calculate Anderson bound of 1D XXZ model ground state energy.
  Use full basis of quantum states without symmetry considerations.
  Use real number elements and symmetric matrix for higher performance.
*/

#ifndef XXZ_1D_ANDERSON_FULL_BASIS_FIXED_NUM_TEST_CPP
#define XXZ_1D_ANDERSON_FULL_BASIS_FIXED_NUM_TEST_CPP

#include <chrono>
#include <cstdio>
#include <stdexcept>

#include "../../hamiltonians_XXZ.hpp"
#include "../include/arlnsmat.h"
#include "../include/arlsmat.h"
#include "../include/arlssym.h"
#include "../matrices/sym/lsmatrxa.h"
#include "../matrices/sym/lsymsol.h"

using std::cout;
using std::endl;
using std::vector;

int main() {
  /*Set parameters sites and Jz.*/
  size_t sites = 20;
  vector<double> Jz;
  for (int i = -8; i < 25; i++) {
    Jz.push_back((double)i / 8.0);
  }
  cout << "Number of sites = " << sites << endl;
  omp_set_num_threads(8);  //Multi-threading when generating matrices
  vector<double> Energy;   //Store the results

  /*Set threads for generating the sparse matrix.
    Use single thread for ARPACK++ since OpenBLAS works better at single thread.*/
  omp_set_num_threads(8);

  for (size_t i = 0; i < Jz.size(); i++) {
    /*Construct polynomial and basis.*/
    SpinHalfPolynomial1D poly = XXZ1D::makeSpinPoly(sites, Jz[i]);
    SpinHalfBasis1D * basis = new SpinHalfBasis1D(sites);
    auto start_basis_init = std::chrono::high_resolution_clock::now();
    basis->init();
    auto end_basis_init = std::chrono::high_resolution_clock::now();
    std::cout << "Basis construction complete!" << std::endl;
    //std::cout << "Basis:" << std::endl << basis->toString() << std::endl;
    //std::cout << "toDecimal:" << (*basis)[0].toDecimal() << (*basis)[1].toDecimal()
    //          << (*basis)[2].toDecimal() << endl;
    const size_t dim = basis->getSize();
    std::cout << "dim = " << dim << std::endl;
    auto duration_basis_init = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_basis_init - start_basis_init);
    cout << "\nInitiating Basis Running Time: " << duration_basis_init.count() << " ms"
         << endl;

    /*Consruct sparse Hamiltonian.*/
    XXZSparseRealHamiltonian * ham = new XXZSparseRealHamiltonian(poly, sites, Jz[i]);
    auto start_matrix_init = std::chrono::high_resolution_clock::now();
    ham->createFullBasisMatrix(*basis);
    auto end_matrix_init = std::chrono::high_resolution_clock::now();
    std::cout << "Hamiltonian construction complete!" << std::endl;
    delete basis;
    //std::cout << "Full Basis:\n" << basis.toString() << std::endl;
    auto duration_matrix_init = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_matrix_init - start_matrix_init);
    cout << "\nInitiating Matrix Running Time: " << duration_matrix_init.count() << " ms"
         << endl;

    int nnz = ham->getNumNonZero();
    std::cout << "nnz = " << nnz << std::endl;
    vector<int> irowVec = ham->getIrow();
    vector<int> pcolVec = ham->getPcol();
    vector<double> valAVec = ham->getNzVal();
    const size_t valASize = valAVec.size();
    delete ham;

    //ARluNonSymMatrix<complex<double>, double> A(dim, nnz, valA, irow, pcol);
    ARluSymMatrix<double> A(dim, nnz, valAVec.data(), irowVec.data(), pcolVec.data());

    // Defining what we need: the 3 lowest eigenvalues of A.
    const size_t solutionNum = 3;
    ARluSymStdEig<double> dprob(solutionNum, A, "SA");
    dprob.ChangeMaxit(1000);

    // Finding eigenvalues and eigenvectors.
    auto start_solve = std::chrono::high_resolution_clock::now();
    dprob.FindEigenvectors();
    auto end_solve = std::chrono::high_resolution_clock::now();

    // Record time.
    auto duration_solve =
        std::chrono::duration_cast<std::chrono::milliseconds>(end_solve - start_solve);
    cout << "\nInitiating Basis Running Time: " << duration_basis_init.count() << " ms"
         << endl;
    cout << "\nInitiating Matrix Running Time: " << duration_matrix_init.count() << " ms"
         << endl;
    cout << "\nSolving Running Time: " << duration_solve.count() << " ms" << endl;

    // Printing solution.
    Solution(A, dprob);

    // Compute and print ground state energy.
    double gs = dprob.Eigenvalue(0);
    for (size_t i = 1; i < solutionNum; i++) {
      if (dprob.Eigenvalue(i) < gs) {
        gs = dprob.Eigenvalue(i);
      }
    }
    gs /= (sites - 2);
    cout << "\nGround State Energy = " << gs << endl;
    Energy.push_back(gs);
  }

  // Print to data file.
  std::string directory = "./_Data/";
  std::string filename =
      directory + "E0vsJz_FullBasis_N_" + std::to_string(sites - 2) + ".txt";
  std::ofstream outFile(filename);
  if (!outFile) {
    throw std::runtime_error("Error: Could not open file " + filename);
  }
  outFile << "Jz"
          << " "
          << "Energy" << std::endl;
  for (size_t i = 0; i < Jz.size(); ++i) {
    outFile << Jz[i] << " " << Energy[i] << std::endl;
  }
  cout << "\nData successfully written in File:\n" << filename << endl << endl;
  outFile.close();

  // Exit
  cout << "EXIT SUCCESS" << endl;
  return EXIT_SUCCESS;
}  // main

#endif  //XXZ_1D_ANDERSON_FULL_BASIS_FIXED_NUM_TEST_CPP
