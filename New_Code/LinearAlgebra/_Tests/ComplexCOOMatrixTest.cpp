/*
  Jiazheng Sun
  Updated: Aug 8, 2024

  Test correctness and performance of ComplexCOOMatrix.
*/

#ifndef LA_COMPLEX_COO_MATRIX_TEST_CPP
#define LA_COMPLEX_COO_MATRIX_TEST_CPP

#include <cstdlib>
#include <iostream>

#include "../printTools.hpp"
#include "../sparseCOO_Tem.hpp"

using std::complex;
using std::cout;
using std::endl;

int main(void) {
  cout << "Complex COO Matrix Test" << endl << endl;

  //Constructor tests
  cout << "-------------------------------------------" << endl;
  cout << "Constructor Tests" << endl;
  ComplexCOOMatrix mat1(15, 15);
  mat1.addData(0, 0, complex<double>(0, 0));
  mat1.addData(1, 1, complex<double>(1, 1));
  mat1.addData(0, 1, complex<double>(0, 1));
  mat1.addData(2, 0, complex<double>(2, 0));
  mat1.addData(0, 2, complex<double>(0, 2));
  mat1.addData(1, 0, complex<double>(1, 0));
  ComplexCOOMatrix mat1cp(mat1);
  cout << "mat1 rows = " << LA::size_tVector_toString(mat1.getRows()) << endl;
  cout << "mat1 cols = " << LA::size_tVector_toString(mat1.getCols()) << endl;
  cout << "mat1 data = " << LA::complexVector_toString(mat1.getAllData()) << endl;
  mat1cp += mat1;
  cout << "\nmat1cp += mat1" << endl;
  cout << "mat1cp rows = " << LA::size_tVector_toString(mat1cp.getRows()) << endl;
  cout << "mat1cp cols = " << LA::size_tVector_toString(mat1cp.getCols()) << endl;
  cout << "mat1cp data = " << LA::complexVector_toString(mat1cp.getAllData()) << endl;
  mat1cp -= mat1;
  cout << "\nmat1cp -= mat1" << endl;
  cout << "mat1cp rows = " << LA::size_tVector_toString(mat1cp.getRows()) << endl;
  cout << "mat1cp cols = " << LA::size_tVector_toString(mat1cp.getCols()) << endl;
  cout << "mat1cp data = " << LA::complexVector_toString(mat1cp.getAllData()) << endl;
  mat1cp *= complex<double>(0, 1);
  cout << "\nmat1cp *= complex<double>(0, 1)" << endl;
  cout << "mat1cp rows = " << LA::size_tVector_toString(mat1cp.getRows()) << endl;
  cout << "mat1cp cols = " << LA::size_tVector_toString(mat1cp.getCols()) << endl;
  cout << "mat1cp data = " << LA::complexVector_toString(mat1cp.getAllData()) << endl;

  //Exit
  cout << "\nEXIT SUCCESS" << endl;
  return EXIT_SUCCESS;
}

#endif  //LA_COMPLEX_COO_MATRIX_TEST_CPP
