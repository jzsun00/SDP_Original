/*
  Jiazheng Sun
  Updated: Jul 28, 2024
*/

#ifndef QM_SPIN_HALF_BASIS_TEST_CPP
#define QM_SPIN_HALF_BASIS_TEST_CPP

#include "../spinStates1D.hpp"

int main(void) {
  SpinHalfBasis1D basis(3);
  basis.init();
  std::cout << "Spin Half Basis:\n" << basis.toString() << std::endl;
  return EXIT_SUCCESS;
}

#endif  //QM_SPIN_HALF_BASIS_TEST_CPP
