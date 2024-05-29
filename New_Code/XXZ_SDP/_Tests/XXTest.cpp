#include "../../HardCore/hardCoreConstraints.hpp"
#include "../../HardCore/hardCoreSubspaces.hpp"
#include "../hamiltonians_XXZ.hpp"

using std::cout;
using std::endl;

int main(void) {
  size_t sites = 91;
  std::string fileName = "XX_N_" + std::to_string(sites) + ".dat";
  printMatrixXX1D(sites, fileName);
  /////////////////////////////////////////////////////
  cout << "\nTests pass!" << endl;
  return EXIT_SUCCESS;
}
