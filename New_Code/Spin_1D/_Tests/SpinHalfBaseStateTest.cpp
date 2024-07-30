/*
  Jiazheng Sun
  Updated: Jul 28, 2024
  
  Test implementations of SpinHalfBaseState1D.
*/

#ifndef QM_SPIN_HALF_BASESTATE_TEST_CPP
#define QM_SPIN_HALF_BASESTATE_TEST_CPP

#include "../spinStates1D.hpp"

int main(void) {
  SpinHalfBaseState1D s0;
  SpinHalfBaseState1D s1(10);
  SpinHalfBaseState1D s2(s1);
  std::cout << "s0 = " << s0.toString() << std::endl;
  std::cout << "s1 = " << s1.toString() << std::endl;
  std::cout << "s2 = " << s2.toString() << std::endl;
  std::cout << "Test pass!" << std::endl;
  return EXIT_SUCCESS;
}

#endif  //QM_SPIN_HALF_BASESTATE_TEST_CPP
