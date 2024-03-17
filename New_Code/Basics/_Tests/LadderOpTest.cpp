#include <cassert>

#include "../operators.hpp"

int main(void) {
  /*Constructors tests: default, (index, creatorF), copy*/

  LadderOp op1(0, true);
  LadderOp op2(1, true);
  LadderOp op3;
  LadderOp op4(1, true);
  LadderOp op5(op1);
  LadderOp op6 = op1;
  assert(op2.getIndex() == 1);
  assert(op2.getCreatorF() == true);
  op3.herm();
  assert(op3.getIndex() == std::numeric_limits<int>::min());
  std::cout << "INT_MIN = " << op3.getIndex() << std::endl;
  std::cout << "INT_MAX = " << std::numeric_limits<int>::max() << std::endl;
  assert(op3.getCreatorF() == true);
  std::cout << "Constructors tests pass!" << std::endl;

  /*Operator overloading tests: =, ==, <, >*/

  assert(op5 == op1);
  assert(op6 == op1);
  assert(&op5 != &op1);
  assert(&op6 != &op1);
  assert(op4 == op2);
  assert(&op4 != &op2);
  assert(op1 < op2);
  assert(op2 > op3);
  try {
    op3.herm();
    op2 > op3;
  }
  catch (std::invalid_argument & e) {
    std::cerr << "Caught exception: " << e.what() << std::endl;
  }
  std::cout << "Operator overload tests pass!" << std::endl;

  /*toString method tests*/

  assert(op1.toString() == "a_{0}{+}");
  assert(op2.toString() == "a_{1}{+}");
  op4.herm();
  assert(op4.toString() == "a_{1}");
  std::cout << "op2 = " << op2.toString() << std::endl;
  std::cout << "op4 = " << op4.toString() << std::endl;
  std::cout << "toString method tests pass!" << std::endl;
  std::cout << "Ladder operator tests run succeccfully!" << std::endl;

  return EXIT_SUCCESS;
}
