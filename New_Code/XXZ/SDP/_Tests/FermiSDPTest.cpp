/*
  Jiazheng Sun
  Updated: Aug 2, 2024
*/

#include "../../../Fermion/fermiConstraints.hpp"
#include "../../../Fermion/fermiSubspaces.hpp"
#include "../../hamiltonians_XXZ.hpp"

using std::complex;
using std::cout;
using std::endl;
using std::pair;
using std::vector;

int main(void) {
  size_t sites = 12;
  double Jz = 0;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly1 =
      makeFermiPoly(-2, sites + 2, Jz);
  cout << "sites = " << sites << "\nJz = " << Jz << endl;
  //cout << "\nHamiltonian =\n" << poly1.toString() << endl;
  poly1.normalize();
  //cout << "\nUse Normal Order:\n"
  //     << "Hamiltonian =\n"
  //     << poly1.toString() << endl;
  //////////////////////////////////////////////////////
  cout << "\nNow construct the spaces" << endl;
  Fermi1DOpSubBasis sub2(0, sites - 1, 2);
  Fermi1DOpSubBasis sub4(20, sites - 21, 4);
  sub2.init();
  sub4.init();
  Fermi1DOpBasis basis;
  basis.addSubspace(sub2);
  basis.addSubspace(sub4);
  //basis.addSubspace(sub6);
  cout << "Operator Basis:" << endl;
  //cout << basis.toString() << endl;
  vector<pair<size_t, size_t> > pairs = FermiFindHermPairs(basis);
  //cout << "\nHermitian Conjugate Pairs:\n" << printHermPairs(pairs) << endl;
  /////////////////////////////////////////////////////
  vector<complex<double> > ham = basis.projPoly(poly1);
  cout << "Hamiltonian Vector:" << endl;
  //cout << complexVector_toString(ham) << endl;
  FermiTransVecToReIm(ham, pairs);
  cout << "\nAfter Transform To Real And Imaginary Parts Of Green's "
          "Functions\nHamiltonian Vector:"
       << endl;
  cout << "Identity Constant:\n" << complex_toString(ham[0]) << endl;
  //cout << complexVector_toString(ham) << endl;
  /////////////////////////////////////////////////////
  cout << "\nNow construct constraint operator set" << endl;
  Fermi1DConsBaseSet base1(0, sites - 1, 1);
  Fermi1DConsBaseSet base2(20, sites - 21, 2);
  //HardCore1DConsBaseSet base3(22, sites - 23, 3);
  base1.init();
  base2.init();
  //base3.init();
  Fermi1DConsSet fullSet;
  fullSet.addBaseSet(base1);
  fullSet.addBaseSet(base2);
  //fullSet.addBaseSet(base3);
  cout << "Constraint operator set:" << endl;
  //cout << fullSet.toString() << endl;
  //std::string fileName = "N_" + std::to_string(sites) + "_Jz_" + ".dat";
  std::string fileNameS = "./XXZ_data/N_" + std::to_string(sites) + ".dat-s";
  cout << "\nNow start writing data files" << endl;
  //printMatrixHardCore1D(fullSet, basis, fileName, ham, pairs);
  printSparseMatrixFermi1D(fullSet, basis, fileNameS, ham, pairs);
  /////////////////////////////////////////////////////
  cout << "\nTests pass!" << endl;
  return EXIT_SUCCESS;
}
