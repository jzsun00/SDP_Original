#include <cassert>

#include "../Basics/fermiStates.hpp"

int main(void) {
  FermiFockState fs1;
  FermiFockState fs2(2);
  std::cout << "fs1 = " << fs1.toString() << std::endl;
  std::cout << "fs2 = " << fs2.toString() << std::endl;
  return EXIT_SUCCESS;
}
