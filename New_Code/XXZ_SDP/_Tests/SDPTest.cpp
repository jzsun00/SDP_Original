#include "../../HardCore/hardCoreConstraints.hpp"
#include "../../HardCore/hardCoreSubspaces.hpp"
#include "../hamiltonians_XXZ.hpp"

using std::cout;
using std::endl;

int main(void) {
  size_t sites = 6;
  double Jz = 0.0;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly1 = makePoly(sites, Jz);
  cout << "sites = " << sites << "\nJz = " << Jz << endl;
  //cout << "poly1 = " << poly1.toString() << endl;
  poly1.normalize();
  //cout << "After normalization\n"
  //     << "poly1 = " << poly1.toString() << endl;
  //////////////////////////////////////////////////////
  cout << "\nNow construct the spaces" << endl;
  HardCore1DOpSubBasis sub1(0, sites - 1, 1);
  HardCore1DOpSubBasis sub2(0, sites - 1, 2);
  sub1.init();
  sub2.init();
  HardCore1DOpBasis basis;
  //basis.addSubspace(sub1);
  basis.addSubspace(sub2);
  cout << "Operator Basis:" << endl;
  //cout << basis.toString() << endl;
  /////////////////////////////////////////////////////
  vector<complex<double> > ham = basis.projPoly(poly1);
  cout << "Hamiltonian Vector:" << endl;
  //cout << complexVector_toString(ham) << endl;
  /////////////////////////////////////////////////////
  cout << "\nNow construct constraint operator set" << endl;
  HardCore1DConsBaseSet base1(0, sites - 1, 1);
  base1.init();
  HardCore1DConsSet fullSet;
  fullSet.addBaseSet(base1);
  cout << "Constraint operator set:" << endl;
  //cout << fullSet.toString() << endl;
  printMatrixHardCore1D(fullSet, basis, "test.dat", ham);
  /////////////////////////////////////////////////////
  cout << "Tests pass!" << endl;
  return EXIT_SUCCESS;
}
