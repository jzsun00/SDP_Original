/*
  Jiazheng Sun
  Updated: Jul 22, 2024

  Test correctness and performance of DoubleDenseMatrix.
*/

#ifndef LA_DOUBLE_DENSE_MATRIX_TEST_CPP
#define LA_DOUBLE_DENSE_MATRIX_TEST_CPP

#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "../dense.hpp"
#include "../dense_Tem.cpp"
#include "../printTools.hpp"

using std::cout;
using std::endl;

int main(void) {
  cout << "Double Dense Matrix Test" << endl << endl;
  std::srand(std::time(0));

  //Constructor tests
  cout << "-------------------------------------------" << endl;
  cout << "Constructor Tests" << endl;
  size_t nrows1 = 3;
  size_t ncols1 = 4;
  size_t nrows2 = 4;
  size_t ncols2 = 2;
  size_t nrows3 = 5;
  size_t ncols3 = 5;
  DoubleDenseMatrix mat1(nrows1, ncols1);
  DoubleDenseMatrix mat2(nrows2, ncols2);
  DoubleDenseMatrix mat3(nrows3, ncols3);
  mat1.fillRandomNum(20);
  mat2.fillRandomNum(20);
  mat3.fillRandomNum(20);
  cout << "nrows1 = " << nrows1 << ",  ncols1 = " << ncols1 << endl;
  cout << "mat1 = " << endl << mat1.toString() << endl << endl;
  cout << "nrows2 = " << nrows2 << ",  ncols2 = " << ncols2 << endl;
  cout << "mat2 = " << endl << mat2.toString() << endl << endl;
  DoubleDenseMatrix mat2copy(mat2);
  cout << "mat2copy = " << endl << mat2copy.toString() << endl << endl;

  //Get information methods tests
  cout << "-------------------------------------------" << endl;
  cout << "Get Information Methods Tests" << endl;
  cout << "nrows1 = " << mat1.getNrows() << ",  "
       << "ncols1 = " << mat1.getNcols() << endl;
  cout << "mat1[0][1] = " << mat1(0, 1) << ",  mat1[1][1] = " << mat1(1, 1)
       << ",  mat1[1][2] = " << mat1(1, 2) << ",  mat1[2][1] = " << mat1(2, 1) << endl
       << endl;

  //Multiplication tests
  cout << "-------------------------------------------" << endl;
  cout << "Multiplication Tests" << endl;
  DoubleDenseMatrix mat12(mat1 * mat2);
  cout << "mat12 = mat1 * mat2 =  " << endl << mat12.toString() << endl << endl;

  //Solve eigenvalue tests
  cout << "-------------------------------------------" << endl;
  cout << "Solve Eigenvalue Tests" << endl;
  cout << "nrows3 = " << nrows3 << ",  ncols3 = " << ncols3 << endl;
  cout << "mat3 = " << endl << mat3.toString() << endl << endl;
  std::vector<double> real(nrows3);
  std::vector<double> imag(nrows3);
  mat3.solveEigen(real, imag);
  cout << "real =" << endl << LA::doubleVector_toString(real) << endl;
  cout << "imag =" << endl << LA::doubleVector_toString(imag) << endl << endl;

  //Exit
  cout << "-------------------------------------------" << endl;
  cout << "\nTests pass!" << endl;
  return EXIT_SUCCESS;
}

#endif  //LA_DOUBLE_DENSE_MATRIX_TEST_CPP
