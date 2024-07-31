#include "../hamiltonians_XXZ.hpp"

int main(void) {
  size_t sites = 3;
  double Jz = 0.2;
  std::cout << "sites = " << sites << ", Jz = " << Jz << std::endl;
  //SpinHalfPolynomial1D poly = makePoly(sites, Jz);
  //std::cout << "poly = " << poly.toString() << std::endl;
  //SpinHalfBasis1D basis(sites);
  //basis.init();
  //XXZSparseHamiltonian ham(poly, sites, Jz);
  //std::cout << "Creating Hamiltonian, now printing the Hamiltonian" << std::endl;
  //std::cout << "poly = " << ham.getPoly().toString() << std::endl;
  //std::cout << "Now constructing the matrix";
  //ham.createMatrix(basis);
  //std::cout << "Printing the Hamiltonian\n" << ham.toString() << std::endl;
  //std::cout << "Test pass!" << std::endl;
  return EXIT_SUCCESS;
}
