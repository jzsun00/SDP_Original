#include "../fermiOperators.hpp"

using std::cout;
using std::endl;

int main(void) {
  Fermi1DLadderOp op1(1, true);
  Fermi1DLadderOp op2(2, false);
  Fermi1DLadderOp op3(3, false);
  Fermi1DLadderOp op4(4, true);
  Fermi1DLadderOp op5(5, true);
  FermiMonomial<Fermi1DLadderOp> mn1(op1);
  FermiMonomial<Fermi1DLadderOp> mn2(op2);
  /*Constructors tests*/
  cout << "Constructors tests" << endl;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly0;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly1(mn1);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly2(complex<double>(-1, 0), mn2);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly1cp(poly1);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly2cp(poly2);
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
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly1cp2 = poly1;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly2cp2 = poly2;
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
  return EXIT_SUCCESS;
}
