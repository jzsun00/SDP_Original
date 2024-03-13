#include "../spinOperators1D.hpp"

int main(void) {
  SpinHalfBaseState b1(5);
  SpinUDHalfOp Op2(2, true);
  SpinHalfState s1 = Op2 * b1;
  std::cout << "s1 = " << s1.toString() << std::endl;
  std::cout << "Test pass!" << std::endl;
  return EXIT_SUCCESS;
}
