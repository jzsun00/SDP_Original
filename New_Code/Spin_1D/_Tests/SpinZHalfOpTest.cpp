#include "../spinOperators1D.hpp"

int main(void) {
  SpinHalfBaseState b1(5);
  SpinZHalfOp Sz1(1);
  SpinHalfState s1 = Sz1 * b1;
  std::cout << "b1 = " << b1.toString() << std::endl;
  std::cout << "Sz1 = " << Sz1.toString() << std::endl;
  std::cout << "s1 = " << s1.toString() << std::endl;
  SpinHalfState s2 = Sz1 * s1;
  std::cout << "s2 = " << s2.toString() << std::endl;
  std::cout << "Test pass!" << std::endl;
  return EXIT_SUCCESS;
}
