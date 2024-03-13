#include "../spinOperators1D.hpp"

int main(void) {
  SpinZHalfOp Sz1(1);
  SpinUDHalfOp Sud1(2, true);
  std::cout << "Sz1 = " << Sz1.toString() << std::endl;
  std::cout << "Sud1 = " << Sud1.toString() << std::endl;
  SpinHalfOp * Sz1cp = new SpinZHalfOp(Sz1);
  SpinHalfMonomial mn1(Sz1cp);
  SpinHalfOp * Sud1cp = new SpinUDHalfOp(Sud1);
  mn1.addOp(Sud1cp);
  std::cout << "mn1 = " << mn1.toString() << std::endl;
  std::cout << "size(mn1) = " << mn1.getSize() << std::endl;
  SpinHalfBaseState base1(5);
  std::cout << "base1 = " << base1.toString() << std::endl;
  SpinHalfState state1 = mn1 * base1;
  std::cout << "mn1 * base1 = " << state1.toString() << std::endl;
  SpinZHalfOp Sz2(2);
  SpinHalfOp * Sz2cp = new SpinZHalfOp(Sz2);
  SpinUDHalfOp Sud0(0, true);
  SpinUDHalfOp * Sud0cp = new SpinUDHalfOp(Sud0);
  SpinHalfMonomial mn2;
  mn2.addOp(Sz2cp);
  mn2.addOp(Sud0cp);
  SpinHalfState state2 = mn2 * state1;
  std::cout << "mn2 = " << mn2.toString() << std::endl;
  std::cout << "size(mn2) = " << mn2.getSize() << std::endl;
  std::cout << "mn2 * mn1 * base1 = " << state2.toString() << std::endl;
  /////////////////////////
  SpinHalfPolynomial poly1;
  std::cout << "Created empty poly1" << std::endl;
  //std::cout << "poly1 = " << poly1.toString() << std::endl;
  poly1 += mn1;
  std::cout << "Added mn1" << std::endl;
  std::cout << "poly1 = " << poly1.toString() << std::endl;
  poly1 += mn2;
  std::cout << "poly1 = mn1 + mn2 = " << poly1.toString() << std::endl;
  SpinHalfState state3 = poly1 * base1;
  std::cout << "poly1 * base1 = " << state3.toString() << std::endl;
  /////////////////////////
  mn1.herm();
  std::cout << "mn1{+} = " << mn1.toString() << std::endl;
  delete Sz1cp;
  delete Sud1cp;
  delete Sz2cp;
  delete Sud0cp;
  return EXIT_SUCCESS;
}
