/*
  Jiazheng Sun
  Updated: Jun 21, 2024

  Calculate Anderson bound of ground state energy of 1D XXZ model.
*/

#include <cstddef>
#include <cstdio>
#include <fstream>

#include "../../hamiltonians_XXZ.hpp"
#include "../include/arcomp.h"
#include "../include/arlnsmat.h"
#include "../include/arlscomp.h"
#include "../matrices/complex/lcompsol.h"

using std::cout;
using std::endl;

int main() {
  /*Set parameters sites and Jz.*/
  size_t sites = 26;
  vector<double> Jz;
  //Jz.push_back(-1);
  for (int i = 1; i < 25; i++) {
    Jz.push_back((double)i / 8.0);
  }
  cout << "Number of sites = " << sites << endl;
  omp_set_num_threads(8);  //Multi-threading when generating matrices
  vector<double> Energy;   //Store the results

  for (size_t i = 0; i < Jz.size(); i++) {
    cout << "Jz = " << Jz[i] << endl << endl;
    /*Construct polynomial and basis.*/
    SpinHalfPolynomial1D poly = makePoly(sites, Jz[i]);
    SpinHalfBasis1D basis(sites);
    auto start_basis_init = std::chrono::high_resolution_clock::now();
    basis.init(0);
    auto end_basis_init = std::chrono::high_resolution_clock::now();
    //std::cout << "Basis construction complete!" << std::endl;
    //std::cout << "Basis:" << std::endl << basis.toString() << std::endl;
    size_t dim = basis.getSize();
    //std::cout << "dim = " << dim << std::endl;

    /*Consruct sparse Hamiltonian.*/
    XXZSparseHamiltonian * ham = new XXZSparseHamiltonian(poly, sites, Jz[i]);
    auto start_matrix_init = std::chrono::high_resolution_clock::now();
    ham->createMatrix(basis);
    auto end_matrix_init = std::chrono::high_resolution_clock::now();
    //std::cout << "Hamiltonian construction complete!" << std::endl;
    //std::cout << "Full Basis:\n" << basis.toString() << std::endl;

    int nnz = ham->getNumNonZero();
    vector<int> irowVec = ham->getIrow();
    vector<int> pcolVec = ham->getPcol();
    vector<complex<double> > valAVec = ham->getNzVal();
    delete ham;
    //cout << "irow = " << intVector_toString(irowVec) << endl;
    //cout << "pcol = " << intVector_toString(pcolVec) << endl;
    //cout << "valA = " << complexVector_toString(valAVec) << endl;

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
    //cout << "\nInitiating Basis Running Time: " << duration_basis_init.count() << " ms"
    //     << endl;
    //cout << "\nInitiating Matrix Running Time: " << duration_matrix_init.count() << " ms"
    //     << endl;
    //cout << "\nSolving Running Time: " << duration_solve.count() << " ms" << endl;

    // Printing solution.
    Solution(A, dprob);
    double gs = dprob.Eigenvalue(0).real();
    for (size_t i = 1; i <= 2; i++) {
      if (dprob.Eigenvalue(i).real() < gs) {
        gs = dprob.Eigenvalue(i).real();
      }
    }
    gs /= (sites - 2);
    cout << "\nGround State Energy = " << gs << endl;
    Energy.push_back(gs);
  }

  // Print to data file.
  std::string directory = "../_Data/";
  std::string filename = "E0vsJz_N_" + std::to_string(sites - 2) + "_pos.txt";
  std::ofstream outFile(filename);
  if (!outFile) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
  }
  outFile << "Jz"
          << "\t"
          << "Energy" << std::endl;
  for (size_t i = 0; i < Jz.size(); ++i) {
    outFile << Jz[i] << "\t" << Energy[i] << std::endl;
  }
  cout << "\nData successfully written in" << filename << endl;
  outFile.close();

}  // main
