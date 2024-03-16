#include "../hamiltonians_XXZ.hpp"

int main(void) {
  size_t sites = 3;
  size_t dim = std::pow(2, sites);
  SpinZHalfOp Sz1(1);
  SpinHalfOp * Sz1cp = new SpinZHalfOp(Sz1);
  SpinHalfMonomial mn1(Sz1cp);
  SpinHalfPolynomial poly1;
  poly1 += mn1;
  std::cout << "poly1 = " << poly1.toString() << std::endl;
  XXZFullHamiltonian ham(poly1, dim);
  SpinHalfBasis basis(sites);
  basis.init();
  std::cout << "Full Basis:\n" << basis.toString() << std::endl;
  ham.createMatrix(basis);
  std::cout << "Hamiltonian Matrix:\n" << ham.toString() << std::endl;
  SpinHalfPolynomial XXZPoly = ham.makePoly(sites, 0.2);
  std::cout << "Polynomial for the XXZ model:\n" << XXZPoly.toString() << std::endl;
  XXZFullHamiltonian ham2(XXZPoly, sites, 0.2);
  ham2.createMatrix(basis);
  std::cout << "Hamiltonian Matrix:\n" << ham2.toString() << std::endl;
  SpinHalfBaseState base0 = basis[0];
  std::cout << "base0 = " << base0.toString() << std::endl;
  SpinHalfState state1 = XXZPoly * base0;
  std::cout << "XXZPoly * base1 = " << state1.toString() << std::endl;
  XXZSparseHamiltonian hamS(XXZPoly, sites, 0.2);
  hamS.createMatrix(basis);
  std::cout << "Sparse Hamiltonian Matrix:\n" << hamS.toString() << std::endl;
  delete Sz1cp;
  return EXIT_SUCCESS;
}
