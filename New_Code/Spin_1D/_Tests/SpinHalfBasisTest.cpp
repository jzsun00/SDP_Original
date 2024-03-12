#include "../spinStates1D.hpp"

int main(void) {
  SpinHalfBasis basis(3);
  basis.init();
  std::cout << "Spin Half Basis:\n" << basis.toString() << std::endl;
  return EXIT_SUCCESS;
}
