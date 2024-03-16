#include <cassert>

#include "../Basics/operators.hpp"

int main(void) {
  LadderOp op0(0, true);
  LadderOp op1(1, true);
  LadderOp op2(2, true);
  /*Constructors tests: default, (LadderOp & Op), copy*/
  Monomial<LadderOp> mono1;
  Monomial<LadderOp> mono2(op0);
  Monomial<LadderOp> mono3(mono2);
  std::cout << "mono1 = " << mono1.toString() << std::endl;
  std::cout << "mono2 = " << mono2.toString() << std::endl;
  std::cout << "mono3 = " << mono3.toString() << std::endl;
  assert(mono1.getSize() == 0);
  assert(mono2.getSize() == 1);
  assert(mono3.getSize() == 1);
  assert(mono2[0] == op0);
  assert(mono3[0] == op0);
  std::cout << "Constructors tests pass!" << std::endl;
  /*Operator overloading tests: =, ==, [], *=*/
  Monomial<LadderOp> mono4 = mono2;
  std::cout << "mono4 = " << mono4.toString() << std::endl;
  assert(mono4 == mono2);
  assert(mono3 == mono2);
  assert(!(mono1 == mono2));
  assert(&mono3 != &mono2);
  assert(&mono4 != &mono2);
  mono4 *= op1;
  Monomial<LadderOp> mono5 = mono4;
  std::cout << "mono4 * op1 = " << mono4.toString() << std::endl;
  std::cout << "mono5 = " << mono5.toString() << std::endl;
  mono4 *= mono5;
  std::cout << "mono4 * mono5 = " << mono4.toString() << std::endl;
  std::cout << "mono4[0] = " << mono4[0].toString() << std::endl;
  std::cout << "mono4[1] = " << mono4[1].toString() << std::endl;
  std::cout << "mono4[2] = " << mono4[2].toString() << std::endl;
  std::cout << "mono4[3] = " << mono4[3].toString() << std::endl;
  mono4 *= op2;
  std::cout << "mono4 * op2 = " << mono4.toString() << std::endl;
  mono4.herm();
  std::cout << "mono4{+} = " << mono4.toString() << std::endl;
  std::cout << "Operator overload tests pass!" << std::endl;
  std::cout << "Ladder operator tests run succeccfully!" << std::endl;
  return EXIT_SUCCESS;
}
