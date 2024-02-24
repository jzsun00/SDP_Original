#include <cassert>

#include "../Basics/operators.hpp"

int main(void) {
  LadderOp op1(0, true);
  LadderOp op2(1, true);
  Monomial mn1(op1);
  Monomial mn2(op2);
  Polynomial pl1(mn1);
  std::cout << "mn1 = " << mn1.toString() << std::endl;
  std::cout << "pl1 = " << pl1.toString() << std::endl;
  std::cout << "mn2 = " << mn2.toString() << std::endl;
  pl1 += pair<complex<double>, Monomial>(complex<double>(1, 0), mn2);
  std::cout << "pl1 + mn2 = " << pl1.toString() << std::endl;
  pl1 += pair<complex<double>, Monomial>(complex<double>(1, 0), mn1);
  std::cout << "pl1 + mn2 + mn1 = " << pl1.toString() << std::endl;
  Polynomial pl2(pl1);
  std::cout << "pl2 = " << pl2.toString() << std::endl;
  pl1 += pl2;
  std::cout << "pl1 + pl2 = " << pl1.toString() << std::endl;
  pl1 *= pair<complex<double>, Monomial>(complex<double>(1, 0), mn1);
  std::cout << "pl1 * mn1 = " << pl1.toString() << std::endl;
  pl1 *= pl2;
  std::cout << "pl1 * pl2 = " << pl1.toString() << std::endl;
  return EXIT_SUCCESS;
}
