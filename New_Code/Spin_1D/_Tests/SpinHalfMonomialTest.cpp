#include "../spinOperators1D.hpp"

int main(void) {
  SpinHalfOp Sz1(1);
  SpinHalfOp Sud1(2, true);
  std::cout << "Sz1 = " << Sz1.toString() << std::endl;
  std::cout << "Sud1 = " << Sud1.toString() << std::endl;
  SpinHalfMonomial mn1(Sz1);
  mn1 *= Sud1;
  std::cout << "mn1 = " << mn1.toString() << std::endl;
  std::cout << "size(mn1) = " << mn1.getSize() << std::endl;
  SpinHalfBaseState base1(5);
  std::cout << "base1 = " << base1.toString() << std::endl;
  SpinHalfState state1 = mn1 * base1;
  std::cout << "mn1 * base1 = " << state1.toString() << std::endl;
  SpinHalfOp Sz2(2);
  SpinHalfOp Sud0(0, true);
  SpinHalfMonomial mn2;
  mn2 *= Sz2;
  mn2 *= Sud0;
  SpinHalfState state2 = mn2 * state1;
  std::cout << "mn2 = " << mn2.toString() << std::endl;
  std::cout << "size(mn2) = " << mn2.getSize() << std::endl;
  std::cout << "mn2 * mn1 * base1 = " << state2.toString() << std::endl;
  mn1.herm();
  std::cout << "mn1{+} = " << mn1.toString() << std::endl;
  return EXIT_SUCCESS;
}
