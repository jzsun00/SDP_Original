/*
  Jiazheng Sun
  Updated: Jun 12, 2024

  Test correctness and performance of DoubleDenseMatrix.
*/

#ifndef DOUBLEDENSEMATRIXTEST_CPP
#define DOUBLEDENSEMATRIXTEST_CPP

#include <cstddef>
#include <cstdlib>

#include "../dense.hpp"
#include "../dense_Tem.cpp"

using std::cout;
using std::endl;

int main(void) {
  cout << "Double Dense Matrix Test" << endl << endl;
  size_t nrows = 2;
  size_t ncols = 2;
  DoubleDenseMatrix mat1(nrows, ncols);
  mat1.setData(0, 0, 1);
  mat1.setData(0, 1, 2);
  mat1.setData(1, 0, 3);
  mat1.setData(1, 1, 4);
  cout << "mat1 = " << endl << mat1.toString() << endl << endl;
  DoubleDenseMatrix mat2(mat1 * mat1);
  cout << "mat2 = mat1 * mat1 =  " << endl << mat2.toString() << endl << endl;
  //Exit
  cout << "Tests pass!" << endl;
  return EXIT_SUCCESS;
}

#endif  //DOUBLEDENSEMATRIXTEST_CPP
