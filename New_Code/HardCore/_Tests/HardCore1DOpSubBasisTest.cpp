#include "../hardCoreSubspaces.hpp"

using std::cout;
using std::endl;

int main(void) {
  HardCore1DOpSubBasis basis1(0, 2, 2);
  basis1.init();
  cout << basis1.toString() << endl;
  HardCore1DLadderOp op0(0, true);
  HardCore1DLadderOp op1(1, true);
  HardCoreMonomial<HardCore1DLadderOp> mn1(op0);
  mn1 *= op1;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly1(mn1);
  cout << "poly1 = " << poly1.toString() << endl;
  vector<complex<double> > proj1 = basis1.projPoly(poly1);
  cout << "Projection1 = " << endl;
  cout << complexVector_toString(proj1) << endl;
  cout << "Tests pass!" << endl;
  return EXIT_SUCCESS;
}
