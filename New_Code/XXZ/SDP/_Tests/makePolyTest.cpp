#include "../../../HardCore/hardCoreSubspaces.hpp"
#include "../../hamiltonians_XXZ.hpp"

using std::complex;
using std::cout;
using std::endl;
using std::vector;

int main(void) {
  size_t sites = 3;
  double Jz = 0.2;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly1 =
      makeHardCorePoly(sites, Jz);
  cout << "sites = " << sites << "\nJz = " << Jz << endl;
  cout << "poly1 = " << poly1.toString() << endl;
  poly1.normalize();
  cout << "After normalization\n"
       << "poly1 = " << poly1.toString() << endl;
  cout << "\nNow construct the spaces" << endl;
  HardCore1DOpSubBasis sub2(0, 3, 2);
  HardCore1DOpSubBasis sub4(0, 3, 4);
  sub2.init(true);
  sub4.init(true);
  cout << sub2.toString() << endl;
  cout << sub4.toString() << endl;
  HardCore1DOpBasis basis;
  basis.addSubspace(sub2);
  basis.addSubspace(sub4);
  cout << "poly1 = " << poly1.toString() << endl;
  cout << "\nCombine two subspaces:" << endl;
  cout << basis.toString() << endl;
  vector<complex<double> > proj1 = basis.projPoly(poly1);
  cout << complexVector_toString(proj1) << endl;
  cout << "Tests pass!" << endl;
  return EXIT_SUCCESS;
}
