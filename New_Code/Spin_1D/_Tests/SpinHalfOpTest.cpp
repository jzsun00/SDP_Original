/*
  Jiazheng Sun
  Updated: Jul 28, 2024
*/

#ifndef QM_SPIN_HALF_OP_TEST_CPP
#define QM_SPIN_HALF_OP_TEST_CPP

#include "../spinOperators1D.hpp"

int main(void) {
  SpinHalfBaseState1D b1(5);
  SpinHalfOp1D Sz1(1);
  SpinHalfState1D s1 = Sz1 * b1;
  std::cout << "b1 = " << b1.toString() << std::endl;
  std::cout << "Sz1 = " << Sz1.toString() << std::endl;
  std::cout << "s1 = " << s1.toString() << std::endl;
  SpinHalfState1D s2 = Sz1 * s1;
  std::cout << "s2 = " << s2.toString() << std::endl;
  std::cout << "Test pass!" << std::endl;
  return EXIT_SUCCESS;
}

#endif  //QM_SPIN_HALF_OP_TEST_CPP
