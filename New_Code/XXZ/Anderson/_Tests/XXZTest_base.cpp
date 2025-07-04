/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE LCompReg.cc.
   Example program that illustrates how to solve a complex standard
   eigenvalue problem in regular mode using the ARluCompStdEig class.

   1) Problem description:

      In this example we try to solve A*x = x*lambda in regular mode,
      where A is obtained from the standard central difference
      discretization of the convection-diffusion operator
                    (Laplacian u) + rho*(du / dx)
      on the unit square [0,1]x[0,1] with zero Dirichlet boundary
      conditions.

   2) Data structure used to represent matrix A:

      {nnz, irow, pcol, A}: matrix A data in CSC format.

   3) Included header files:

      File             Contents
      -----------      ---------------------------------------------
      lcmatrxa.h       CompMatrixA, a function that generates matrix
                       A in CSC format.
      arlnsmat.h       The ARluNonSymMatrix class definition.
      arlscomp.h       The ARluCompStdEig class definition.
      lcompsol.h       The Solution function.
      arcomp.h         The "arcomplex" (complex) type definition.

   4) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include <cstdio>
#include <string>
#include <vector>

#include "../../../XXZ/hamiltonians_XXZ.hpp"
#include "arcomp.h"
#include "arlnsmat.h"
#include "arlscomp.h"
#include "lcompsol.h"

using std::complex;
using std::vector;

int main() {
  // Defining variables;

  //int nx;
  //int n;                     // Dimension of the problem.
  //int nnz;                   // Number of nonzero elements in A.
  //int * irow;                // pointer to an array that stores the row
  // indices of the nonzeros in A.
  //int * pcol;                // pointer to an array of pointers to the
  // beginning of each column of A in valA.
  //arcomplex<double> * valA;  // pointer to an array that stores the
  // nonzero elements of A.

  // Creating a complex matrix.

  //nx = 10;
  //n = nx * nx;
  //CompMatrixA(nx, nnz, valA, irow, pcol);
  //ARluNonSymMatrix<arcomplex<double>, double> A(n, nnz, valA, irow, pcol);

  size_t sites = 3;
  double Jz = 0.2;
  size_t dim = std::pow(2, sites);
  SpinHalfBasis1D basis(sites);
  //std::cout << "Full Basis:\n" << basis.toString() << std::endl;
  //XXZFullHamiltonian ham0;
  //SpinHalfPolynomial XXZPoly = ham0.makePoly(sites, Jz);

  int n = 8;
  int nnz = 12;
  int * irow = new int[nnz];
  int * pcol = new int[n + 1];
  complex<double> * valA = new complex<double>[nnz];
  irow[0] = 0;
  irow[1] = 2;
  irow[2] = 1;
  irow[3] = 2;
  irow[4] = 4;
  irow[5] = 5;
  irow[6] = 2;
  irow[7] = 3;
  irow[8] = 5;
  irow[9] = 6;
  irow[10] = 5;
  irow[11] = 7;
  pcol[0] = 0;
  pcol[1] = 1;
  pcol[2] = 2;
  pcol[3] = 5;
  pcol[4] = 6;
  pcol[5] = 7;
  pcol[6] = 10;
  pcol[7] = 11;
  pcol[8] = 12;
  valA[0] = complex<double>(0.1, 0.0);
  valA[1] = complex<double>(0.5, 0.0);
  valA[2] = complex<double>(0.5, 0.0);
  valA[3] = complex<double>(-0.1, 0.0);
  valA[4] = complex<double>(0.5, 0.0);
  valA[5] = complex<double>(0.5, 0.0);
  valA[6] = complex<double>(0.5, 0.0);
  valA[7] = complex<double>(0.5, 0.0);
  valA[8] = complex<double>(-0.1, 0.0);
  valA[9] = complex<double>(0.5, 0.0);
  valA[10] = complex<double>(0.5, 0.0);
  valA[11] = complex<double>(0.1, 0.0);

  ARluNonSymMatrix<complex<double>, double> A(n, nnz, valA, irow, pcol);

  // Defining what we need: the four eigenvectors of A with largest magnitude.

  ARluCompStdEig<double> dprob(4L, A);

  // Finding eigenvalues and eigenvectors.

  dprob.FindEigenvectors();

  // Printing solution.

  Solution(A, dprob);

}  // main
