#include "../../HardCore/hardCoreConstraints.hpp"
#include "../../HardCore/hardCoreSubspaces.hpp"
#include "../hamiltonians_XXZ.hpp"

using std::cout;
using std::endl;

int main(void) {
  size_t sites = 9;
  double Jz = 0;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly1 = makePoly(sites, Jz);
  cout << "sites = " << sites << "\nJz = " << Jz << endl;
  cout << "\nHamiltonian =\n" << poly1.toString() << endl;
  poly1.normalize();
  cout << "\nUse Normal Order:\n"
       << "Hamiltonian =\n"
       << poly1.toString() << endl;
  //////////////////////////////////////////////////////
  cout << "\nNow construct the spaces" << endl;
  HardCore1DOpSubBasis sub1(0, sites - 1, 1);
  HardCore1DOpSubBasis sub2(0, sites - 1, 2);
  HardCore1DOpSubBasis sub4(0, sites - 1, 4);
  sub1.init();
  sub2.init();
  //sub4.init();
  HardCore1DOpBasis basis;
  //basis.addSubspace(sub1);
  basis.addSubspace(sub2);
  //basis.addSubspace(sub4);
  cout << "Operator Basis:" << endl;
  cout << basis.toString() << endl;
  vector<pair<size_t, size_t> > pairs = findHermPairs(basis);
  cout << "\nHermitian Conjugate Pairs:\n" << printHermPairs(pairs) << endl;
  /////////////////////////////////////////////////////
  vector<complex<double> > ham = basis.projPoly(poly1);
  cout << "Hamiltonian Vector:" << endl;
  cout << complexVector_toString(ham) << endl;
  transVecToReIm(ham, pairs);
  cout << "\nAfter Transform To Real And Imaginary Parts Of Green's "
          "Functions\nHamiltonian Vector:"
       << endl;
  cout << complexVector_toString(ham) << endl;
  /////////////////////////////////////////////////////
  cout << "\nNow construct constraint operator set" << endl;
  HardCore1DConsBaseSet base1(0, sites - 1, 1);
  HardCore1DConsBaseSet base2(0, sites - 1, 2);
  base1.init();
  //base2.init();
  HardCore1DConsSet fullSet;
  fullSet.addBaseSet(base1);
  //fullSet.addBaseSet(base2);
  cout << "Constraint operator set:" << endl;
  cout << fullSet.toString() << endl;
  std::string fileName = "N_" + std::to_string(sites) + ".dat";
  printMatrixHardCore1D(fullSet, basis, fileName, ham, pairs);
  /////////////////////////////////////////////////////
  cout << "\nTests pass!" << endl;
  return EXIT_SUCCESS;
}
