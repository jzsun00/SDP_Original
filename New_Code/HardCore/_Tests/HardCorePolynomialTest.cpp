#include "../hardCoreOperators.hpp"

using std::cout;
using std::endl;

int main(void) {
  HardCore1DLadderOp op1(1, true);
  HardCore1DLadderOp op2(2, false);
  HardCore1DLadderOp op3(3, false);
  HardCore1DLadderOp op4(4, true);
  HardCore1DLadderOp op5(5, true);
  HardCoreMonomial<HardCore1DLadderOp> mn1(op1);
  HardCoreMonomial<HardCore1DLadderOp> mn2(op2);
  /*Constructors tests*/
  cout << "Constructors tests" << endl;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly0;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly1(mn1);
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly2(complex<double>(-1, 0),
                                                                  mn2);
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly1cp(poly1);
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly2cp(poly2);
  cout << "Default constructor: "
       << "poly0 = " << poly0.toString() << endl;
  cout << "Explicit constructor: "
       << "poly1 = " << poly1.toString() << endl;
  cout << "Explicit constructor: "
       << "poly2 = " << poly2.toString() << endl;
  cout << "Copy constructor: "
       << "poly1cp = " << poly1.toString() << endl;
  cout << "Copy constructor: "
       << "poly2cp = " << poly2.toString() << endl;

  /*Overloaded operators tests*/
  cout << "\nOverloaded operators tests" << endl;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly1cp2 = poly1;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly2cp2 = poly2;
  cout << "Copy operator: "
       << "poly1cp2 = " << poly1.toString() << endl;
  cout << "Copy operator: "
       << "poly2cp2 = " << poly2.toString() << endl;
  poly1cp2 *= op5;
  poly2cp2 *= op4;
  cout << "*= operator: (" << poly1.toString() << " *= " << op5.toString()
       << ") = " << poly1cp2.toString() << endl;
  cout << "*= operator: (" << poly2.toString() << " *= " << op4.toString()
       << ") = " << poly2cp2.toString() << endl;
  poly1cp2 *= poly2cp2;
  poly1cp2.eraseZeros();
  cout << "*= operator: (" << poly1cp2.toString() << " *= " << poly2cp2.toString()
       << ") = " << poly1cp2.toString() << endl;

  /*Normalization tests.*/
  cout << "\nNormalization tests." << endl;
  cout << "poly1cp2 = " << poly1cp2.toString() << endl;
  poly1cp2.normalize();
  poly1cp2.eraseNonNorm();
  cout << "poly1cp2 => " << poly1cp2.toString() << endl;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly11(mn1);
  mn1.herm();
  poly11 *= mn1;
  cout << "poly11 = " << poly11.toString() << endl;
  poly11.normalize();
  poly11.eraseNonNorm();
  cout << "poly11 => " << poly11.toString() << endl;
  poly1cp2 *= poly11;
  poly1cp2.eraseZeros();
  cout << "poly1cp2 *= poly11 = " << poly1cp2.toString() << endl;
  poly1cp2.normalize();
  //poly1cp2.eraseNonNorm();
  cout << "poly1cp2 => " << poly1cp2.toString() << endl;
  return EXIT_SUCCESS;
}
