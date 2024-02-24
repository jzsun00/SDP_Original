#include "../Basics/operators.hpp"

#include <cassert>

int LadderOpTest(void) {
  LadderOp * op1 = new LadderOp(0, true);
  LadderOp * op2 = new LadderOp(1, true);
  LadderOp * op3 = new LadderOp();
  LadderOp * op4 = new LadderOp(1, true);
  LadderOp * op5 = new LadderOp(*op1);
  LadderOp * op6 = new LadderOp();
  *op6 = *op1;
  assert(op2->getIndex() == 1);
  assert(op2->getCreatorF() == true);
  op3->herm();
  assert(op3->getIndex() == std::numeric_limits<int>::min());
  assert(op3->getCreatorF() == true);
  assert(*op5 == *op1);
  assert(op5 != op1);
  std::cout << "Constructor tests pass!" << std::endl;
  assert((*op2 == *op4) == true);
  assert((*op2 > *op1) == true);
  assert((*op1 < *op2) == true);
  assert((*op1 == *op2) == false);
  assert(*op6 == *op1);
  assert(op6 != op1);
  try {
    op3->herm();
    bool mid = *op2 > *op3;
    std::cout << mid << std::endl;
  }
  catch (std::invalid_argument & e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
  std::cout << "Operator overload tests pass!" << std::endl;
  op2->herm();
  std::cout << op1->toString() << std::endl;
  std::cout << op2->toString() << std::endl;
  std::cout << "toString method tests pass!" << std::endl;
  delete op1;
  delete op2;
  delete op3;
  delete op4;
  delete op5;
  delete op6;
  std::cout << "The program runs succeccfully!" << std::endl;
  return EXIT_SUCCESS;
}

int MonomialTest(void) {
  LadderOp * op1 = new LadderOp(0, true);
  LadderOp * op2 = new LadderOp(1, true);
  LadderOp * op3 = new LadderOp(2, false);
  Monomial * mn1 = new Monomial();
  Monomial * mn2 = new Monomial(*op1);
  assert(mn1->getSize() == 0);
  assert(mn2->getSize() == 1);
  *mn1 *= *op1;
  *mn1 *= *op2;
  *mn1 *= *op3;
  assert(mn1->getSize() == 3);
  Monomial * mn3 = new Monomial(*mn1);
  Monomial * mn4 = new Monomial();
  *mn4 = *mn3;
  assert(*mn3 == *mn1);
  assert(*mn4 == *mn3);
  assert(mn3 != mn1);
  assert(mn4 != mn3);
  *mn4 *= *mn2;
  std::cout << "mn1 = " << mn1->toString() << std::endl;
  std::cout << "mn2 = " << mn2->toString() << std::endl;
  std::cout << "mn3 = " << mn3->toString() << std::endl;
  std::cout << "mn4 = " << mn4->toString() << std::endl;
  assert((*mn4)[3] == *op1);
  assert((*mn4)[1] == *op2);
  mn2->herm();
  mn4->herm();
  std::cout << "mn2+ = " << mn2->toString() << std::endl;
  std::cout << "mn4+ = " << mn4->toString() << std::endl;
  return EXIT_SUCCESS;
}

int PolynomialTest(void) {
  LadderOp op1(0, true);
  Monomial mn1(op1);
  Polynomial pl1(mn1);
  std::cout << "pl1 = " << pl1.toString() << std::endl;
  return EXIT_SUCCESS;
}

int main(void) {
  assert(LadderOpTest() == EXIT_SUCCESS);
  assert(MonomialTest() == EXIT_SUCCESS);
  assert(PolynomialTest() == EXIT_SUCCESS);
  return EXIT_SUCCESS;
}
