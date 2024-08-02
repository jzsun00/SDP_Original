/*
  Jiazheng Sun
  Updated: Aug 2, 2024
  
  Calculate Anderson bound of ground state energy of 1D XXZ model.
*/

#include <cstddef>
#include <cstdio>
#include <fstream>
#include <string>

#include "../../hamiltonians_XXZ.hpp"
#include "../include/arlnsmat.h"
#include "../include/arlsmat.h"
#include "../include/arlssym.h"
#include "../matrices/sym/lsmatrxa.h"
#include "../matrices/sym/lsymsol.h"
#include "omp.h"

using std::complex;
using std::cout;
using std::endl;
using std::vector;

int main() {
  /*Set parameters sites and Jz.*/
  int sites = 23;
  double Jz = 0.5;
  vector<int> Sztotal;  //Actual Sz^total should divide by two
  cout << "Number of sites = " << sites << endl;
  cout << "Jz = " << Jz << endl;
  omp_set_num_threads(12);  //Multi-threading when generating matrices
  vector<double> Energy;    //Store the results
  for (int i = -sites + 2; i <= sites - 2; i += 2) {
    Sztotal.push_back(i);
  }

  for (size_t i = 0; i < Sztotal.size(); i++) {
    cout << "Sztotal = " << Sztotal[i] << endl << endl;
    /*Construct polynomial and basis.*/
    SpinHalfPolynomial1D poly = makePoly(sites, Jz);
    SpinHalfBasis1D * basis = new SpinHalfBasis1D(sites);
    auto start_basis_init = std::chrono::high_resolution_clock::now();
    basis->init(Sztotal[i]);
    auto end_basis_init = std::chrono::high_resolution_clock::now();
    //std::cout << "Basis construction complete!" << std::endl;
    //std::cout << "Basis:" << std::endl << basis.toString() << std::endl;
    size_t dim = basis->getSize();
    std::cout << "dim = " << dim << std::endl;

    /*Consruct sparse Hamiltonian.*/
    XXZSparseRealHamiltonian * ham = new XXZSparseRealHamiltonian(poly, sites, Jz);
    auto start_matrix_init = std::chrono::high_resolution_clock::now();
    ham->createMatrix(*basis);
    auto end_matrix_init = std::chrono::high_resolution_clock::now();
    std::cout << "Hamiltonian construction complete!" << std::endl;
    delete basis;
    //std::cout << "Full Basis:\n" << basis.toString() << std::endl;

    int nnz = ham->getNumNonZero();
    vector<int> irowVec = ham->getIrow();
    vector<int> pcolVec = ham->getPcol();
    vector<double> valAVec = ham->getNzVal();
    delete ham;
    //cout << "irow = " << intVector_toString(irowVec) << endl;
    //cout << "pcol = " << intVector_toString(pcolVec) << endl;
    //cout << "valA = " << complexVector_toString(valAVec) << endl;

    ARluSymMatrix<double> A(dim, nnz, valAVec.data(), irowVec.data(), pcolVec.data());

    // Defining what we need: the 3 lowest eigenvalues of A.
    ARluSymStdEig<double> dprob(3, A, "SA");
    dprob.ChangeMaxit(1000);

    // Finding eigenvalues and eigenvectors.
    auto start_solve = std::chrono::high_resolution_clock::now();
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
    double gs = dprob.Eigenvalue(0);
    for (size_t i = 1; i < 3; i++) {
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
      directory + "E0vsSztotal_N_" + std::to_string(sites - 2) + "_Jz_0_5" + ".txt";
  std::ofstream outFile(filename);
  if (!outFile) {
    throw std::runtime_error("Error: Could not open file " + filename);
  }
  outFile << "Sztotal"
          << " "
          << "Energy" << std::endl;
  for (size_t i = 0; i < Sztotal.size(); ++i) {
    outFile << Sztotal[i] << " " << Energy[i] << std::endl;
  }
  cout << "\nData successfully written in File:\n" << filename << endl << endl;
  outFile.close();

  // Exit
  cout << "EXIT SUCCESS" << endl;
  return EXIT_SUCCESS;
}  // main
