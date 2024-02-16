#include "operators.hpp"

#include <cassert>

int main(void) {
  LadderOp * op1 = new LadderOp(0, true);
  LadderOp * op2 = new LadderOp(1, true);
  LadderOp * op3 = new LadderOp(2, false);
  assert(op1->getIndex() == 0);
  assert(op1->getCreatorF() == true);
  assert((*op1 < *op2) == true);
  std::cout << "The program runs succeccfully.";
  return EXIT_SUCCESS;
}
