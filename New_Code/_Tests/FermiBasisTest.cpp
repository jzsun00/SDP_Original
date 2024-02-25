#include <cassert>

#include "../Basics/fermiStates.hpp"

int main(void) {
  FermiBasis fb1(3);
  fb1.init();
  std::cout << "fb1:\n" << fb1.toString() << std::endl;
  return EXIT_SUCCESS;
}
