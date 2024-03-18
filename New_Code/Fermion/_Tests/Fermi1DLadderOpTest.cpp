#include <cassert>

#include "../fermiOperators.hpp"

int main(void) {
  Fermi1DLadderOp lad0(1, true);
  std::cout << "lad0 = " << lad0.toString() << std::endl;
  std::cout << "Tests pass!" << std::endl;
  return EXIT_SUCCESS;
}
