#include <cassert>

#include "../Basics/fermiOperators.hpp"

int main(void) {
  size_t N = 15;
  long unsigned int dim = std::pow(2, N);
  std::cout << "dim = " << dim << std::endl;
  FermiBasis fb1(N);
  fb1.init();
  std::cout << "Initialization finished" << std::endl;
  //std::cout << "fb1:\n" << fb1.toString() << std::endl;
  FermiLadderOp flo1(1, true);
  FermiState fs0 = flo1 * fb1[0];
  FermiState fs5 = flo1 * fb1[5];
  //std::cout << flo1.toString() << " * " << fb1[0].toString() << " = " << fs0.toString()
  //<< std::endl;
  //std::cout << flo1.toString() << " * " << fb1[5].toString() << " = " << fs5.toString()
  //<< std::endl;
  vector<vector<complex<double> > > matrix(dim, vector<complex<double> >(dim));

  for (size_t i = 0; i < std::pow(2, N); i++) {
    for (size_t j = 0; j < std::pow(2, N); j++) {
      //std::cout << "i = " << i << ", "
      //      << "j = " << j << std::endl;
      FermiState fs = flo1 * fb1[j];
      matrix[i][j] = innerProduct(fb1[i], fs);
    }
  }
  std::cout << "Matirx filling finished" << std::endl;
  return EXIT_SUCCESS;
}
