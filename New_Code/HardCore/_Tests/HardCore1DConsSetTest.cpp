#include "../hardCoreConstraints.hpp"

using std::cout;
using std::endl;

int main(void) {
  HardCore1DConsBaseSet base1(0, 1, 1);
  HardCore1DConsBaseSet base2(0, 1, 2);
  base1.init();
  base2.init();
  //cout << base1.toString() << endl;
  //cout << base2.toString() << endl;
  HardCore1DConsSet fullSet;
  fullSet.addBaseSet(base1);
  //fullSet.addBaseSet(base2);
  cout << "Constraint operator set:" << endl;
  cout << fullSet.toString() << endl;
  cout << fullSet.getIJPoly(1, 2).toString() << endl;
  HardCore1DOpSubBasis subBasis1(0, 1, 1);
  subBasis1.init();
  HardCore1DOpSubBasis subBasis2(0, 1, 2);
  subBasis2.init();
  HardCore1DOpBasis basis;
  basis.addSubspace(subBasis1);
  basis.addSubspace(subBasis2);
  cout << "Basis:" << endl;
  cout << basis.toString() << endl;
  printMatrixHardCore1D(fullSet, basis);
  cout << "Tests pass!" << endl;
}
