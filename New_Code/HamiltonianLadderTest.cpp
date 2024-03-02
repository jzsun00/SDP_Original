#include <cassert>

#include "../Basics/hamiltonians.hpp"

int main(void) {
  size_t N = 3;
  unsigned long dim = std::pow(2, N);
  FermiLadderOp op(1, true);
  FermiBasis basis(dim);
  basis.init();
  Hamiltonian<FermiLadderOp, FermiBasis> H(op, dim);
  return EXIT_SUCCESS;
}
