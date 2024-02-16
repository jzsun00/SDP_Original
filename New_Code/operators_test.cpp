#include "operators.hpp"

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
    *op2 > *op3;
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

int main() {
  assert(LadderOpTest() == EXIT_SUCCESS);
  return EXIT_SUCCESS;
}
