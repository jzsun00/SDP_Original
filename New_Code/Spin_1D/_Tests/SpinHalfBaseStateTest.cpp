#include "../spinStates1D.hpp"

int main(void) {
  SpinHalfBaseState s0;
  SpinHalfBaseState s1(10);
  SpinHalfBaseState s2(s1);
  std::cout << "s0 = " << s0.toString() << std::endl;
  std::cout << "s1 = " << s1.toString() << std::endl;
  std::cout << "s2 = " << s2.toString() << std::endl;
  std::cout << "Test pass!" << std::endl;
  return EXIT_SUCCESS;
}
