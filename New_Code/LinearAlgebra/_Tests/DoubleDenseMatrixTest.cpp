/*
  Jiazheng Sun
  Updated: Jun 14, 2024

  Test correctness and performance of DoubleDenseMatrix.
*/

#ifndef DOUBLE_DENSE_MATRIX_TEST_CPP
#define DOUBLE_DENSE_MATRIX_TEST_CPP

#include <cstddef>
#include <cstdlib>
#include <ctime>

#include "../dense.hpp"
#include "../dense_Tem.cpp"

using std::cout;
using std::endl;

int main(void) {
  cout << "Double Dense Matrix Test" << endl << endl;
  std::srand(std::time(0));
  //Constructor tests
  cout << "Constructor tests" << endl;
  size_t nrows1 = 3;
  size_t ncols1 = 4;
  size_t nrows2 = 4;
  size_t ncols2 = 2;
  DoubleDenseMatrix mat1(nrows1, ncols1);
  DoubleDenseMatrix mat2(nrows2, ncols2);
  mat1.fillRandomNum();
  mat2.fillRandomNum();
  cout << "nrows1 = " << nrows1 << ",  ncols1 = " << ncols1 << endl;
  cout << "mat1 = " << endl << mat1.toString() << endl << endl;
  cout << "nrows2 = " << nrows2 << ",  ncols2 = " << ncols2 << endl;
  cout << "mat2 = " << endl << mat2.toString() << endl << endl;
  DoubleDenseMatrix mat2copy(mat2);
  cout << "mat2copy = " << endl << mat2copy.toString() << endl << endl;
  //Get information methods tests
  cout << "Get information methods tests" << endl;
  cout << "nrows1 = " << mat1.getNrows() << ",  "
       << "ncols1 = " << mat1.getNcols() << endl;
  cout << "mat1[1][1] = " << mat1(1, 1) << ",  mat1[1][2] = " << mat1(1, 2)
       << ",  mat1[2][1] = " << mat1(2, 1) << endl
       << endl;
  //Multiplication tests
  cout << "Multiplication tests" << endl;
  DoubleDenseMatrix mat12(mat1 * mat2);
  cout << "mat12 = mat1 * mat2 =  " << endl << mat12.toString() << endl << endl;
  //Exit
  cout << "Tests pass!" << endl;
  return EXIT_SUCCESS;
}

#endif  //DOUBLE_DENSE_MATRIX_TEST_CPP
