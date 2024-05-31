#include "../../HardCore/hardCoreConstraints.hpp"
#include "../../HardCore/hardCoreSubspaces.hpp"
#include "../hamiltonians_XXZ.hpp"

using std::cout;
using std::endl;

int main(void) {
  size_t sites = 151;
  std::string fileName = "XX_N_" + std::to_string(sites) + ".dat";
  std::string fileNameS = "XX_N_" + std::to_string(sites) + ".dat-s";
  //printMatrixXX1D(sites, fileName);
  printSparseMatrixXX1D(sites, fileNameS);
  /////////////////////////////////////////////////////
  cout << "\nTests pass!" << endl;
  return EXIT_SUCCESS;
}
