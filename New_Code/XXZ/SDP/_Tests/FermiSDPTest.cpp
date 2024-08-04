/*
  Jiazheng Sun
  Updated: Aug 3, 2024
  
  Use spinless Fermion formalism to compute ground state energy of 1D XXZ Model.
*/

#ifndef XXZ_1D_SDP_FERMI_TEST_CPP
#define XXZ_1D_SDP_FERMI_TEST_CPP

#include "../../../Fermion/fermiConstraints.hpp"
#include "../../../Fermion/fermiSubspaces.hpp"
#include "../../hamiltonians_XXZ.hpp"

using std::complex;
using std::cout;
using std::endl;
using std::pair;
using std::vector;

int main(void) {
  cout << "\n1D XXZ Model Test: SDP Method" << endl << endl;

  /*Set number of sites and Jz.*/
  size_t sites1 = 20;  //First order
  size_t sites2 = 0;   //Second order
  double Jz = 0;
  cout << "sites1 = " << sites1 << "\nsites2 = " << sites2 << "\nJz = " << Jz << endl;

  /*Construct Hamiltonian and convert to normal order.*/
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly1 =
      makeFermiPolyPBC(0, sites1 - 1, Jz);
  poly1.normalize();
  //cout << "H_XXZ = " << poly1.toString() << endl;

  /*Construct operator basis.*/
  cout << "\nNow construct operator basis" << endl;
  Fermi1DOpBasis basis;
  Fermi1DOpSubBasis sub2(0, sites1 - 1, 2);
  sub2.init();
  basis.addSubspace(sub2);
  if (sites2 != 0) {
    Fermi1DOpSubBasis sub4((sites1 - sites2) / 2, (sites1 + sites2) / 2, 4);
    sub4.init();
    basis.addSubspace(sub4);
  }
  cout << "Operator Basis Construction Finished" << endl;
  //cout << "Operator Basis:" << basis.toString() << endl;
  vector<pair<size_t, size_t> > pairs = FermiFindHermPairs(basis);
  //cout << "\nHermitian Conjugate Pairs:\n" << printHermPairs(pairs) << endl;

  /*Compute cost function vector.*/
  cout << "\nNow compute the cost function vector" << endl;
  vector<complex<double> > ham = basis.projPoly(poly1);
  //cout << "Hamiltonian Vector:" << endl << complexVector_toString(ham) << endl;
  FermiTransVecToReIm(ham, pairs);
  //cout << "\nAfter Transform To Real And Imaginary Parts Of Green's "
  //       "Functions\nHamiltonian Vector:"
  //     << endl;
  cout << "Identity Constant:\n" << complex_toString(ham[0]) << endl;
  //cout << complexVector_toString(ham) << endl;

  /*Compute constraint matrices.*/
  cout << "\nNow construct constraint operator set" << endl;
  Fermi1DConsSet fullSet;
  Fermi1DConsBaseSet base1(0, sites1 - 1, 1);
  base1.init();
  fullSet.addBaseSet(base1);
  if (sites2 != 0) {
    Fermi1DConsBaseSet base2((sites1 - sites2) / 2, (sites1 + sites2) / 2, 2);
    base2.init();
    fullSet.addBaseSet(base2);
  }
  cout << "Constraint Set Construction Finished" << endl;
  //cout << "Constraint operator set:\n" << fullSet.toString() << endl;

  /*Print cost function and constraint matrices into file.*/
  std::string fileNameS =
      "./XXZ_data/N_" + std::to_string(sites1) + "_" + std::to_string(sites2) + ".dat-s";
  cout << "\nNow start writing data files" << endl;
  printSparseMatrixFermi1D(fullSet, basis, fileNameS, ham, pairs);
  cout << "\nData successfully written in File:\n" << fileNameS << endl << endl;

  /*Exit*/
  cout << "\nEXIT SUCCESS" << endl;
  return EXIT_SUCCESS;
}

#endif  //XXZ_1D_SDP_FERMI_TEST_CPP
