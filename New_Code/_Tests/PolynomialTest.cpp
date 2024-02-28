#include <cassert>

#include "../Basics/operators.hpp"

int main(void) {
  LadderOp op0(0, true);
  LadderOp op1(1, true);
  LadderOp op2(2, true);
  Monomial mn0(op0);
  Monomial mn1(op1);
  Monomial mn2(op2);
  /*Constructors tests: default, (Monomial const & mn), copy*/
  Polynomial pl0;
  Polynomial pl1(mn1);
  Polynomial pl2(pl1);
  assert(pl0.getSize() == 0);
  assert(pl1.getSize() == 1);
  assert(pl2.getSize() == 1);
  std::cout << "mn1 = " << mn1.toString() << std::endl;
  std::cout << "pl1 = " << pl1.toString() << std::endl;
  std::cout << "pl2 = " << pl2.toString() << std::endl;
  std::cout << "Constructors tests pass!" << std::endl;
  /*Operator overloading tests: =, ==, [], +=, -+, *=*/
  Polynomial pl3 = pl2;
  assert(pl1 == pl2);
  assert(pl2 == pl3);
  assert(&pl1 != &pl2);
  assert(&pl2 != &pl3);
  std::cout << "pl3 = " << pl3.toString() << std::endl;
  pl1 += mn0;
  pl1 += mn2;
  std::cout << "pl1 + mn0 + mn2 = " << pl1.toString() << std::endl;
  pl2 += mn0;
  std::cout << "pl2 + mn0 = " << pl2.toString() << std::endl;
  pl1 -= pl2;
  std::cout << "pl1 - pl2 = " << pl1.toString() << std::endl;
  std::cout << "len(pl1) = " << pl1.getSize() << std::endl;
  pl1.eraseZeros();
  std::cout << "After erasing zeros:" << std::endl;
  std::cout << "pl1 = " << pl1.toString() << std::endl;
  std::cout << "len(p1) = " << pl1.getSize() << std::endl;
  Polynomial pl4 = pl2;
  std::cout << "pl2 = " << pl2.toString() << std::endl;
  std::cout << "pl4 = " << pl4.toString() << std::endl;
  pl4 *= pl2;
  pl4.eraseZeros();
  std::cout << "pl4 * pl2 = " << pl4.toString() << std::endl;
  pl4.herm();
  std::cout << "pl4{+} = " << pl4.toString() << std::endl;
  return EXIT_SUCCESS;
}
