/*
  Jiazheng Sun
  Updated: Aug 6, 2024
  
  Test implementations of FermiPolynomial<MonomialType>.
*/

#ifndef QM_FERMI_POLYNOMIAL_TEST_CPP
#define QM_FERMI_POLYNOMIAL_TEST_CPP

#include "../fermiOperators_Tem.hpp"

using std::complex;
using std::cout;
using std::endl;

int main(void) {
  cout << "Fermi Polynomial Test\n" << endl;

  /*Constructors tests*/
  cout << "Constructors Tests" << endl;
  Fermi1DLadderOp op1(1, true);
  Fermi1DLadderOp op2(2, false);
  Fermi1DLadderOp op3(3, false);
  Fermi1DLadderOp op4(4, true);
  Fermi1DLadderOp op5(5, true);
  FermiMonomial<Fermi1DLadderOp> mn1(op1);
  FermiMonomial<Fermi1DLadderOp> mn2(op2);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly0;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly1(mn1);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly2(complex<double>(-1, 0), mn2);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly1cp(poly1);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly2cp(poly2);
  cout << "Default constructor: "
       << "poly0 = " << poly0.toString() << endl;
  cout << "Explicit constructor(1, 1, true): "
       << "poly1 = " << poly1.toString() << endl;
  cout << "Explicit constructor(-1, 2, false): "
       << "poly2 = " << poly2.toString() << endl;
  cout << "Copy constructor(poly1): "
       << "poly1cp = " << poly1cp.toString() << endl;
  cout << "Copy constructor(poly2): "
       << "poly2cp = " << poly2cp.toString() << endl;

  /*Overloaded operators tests*/
  cout << "\nOverloaded Operators Tests" << endl;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly1cp2 = poly1;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly2cp2 = poly2;
  cout << "Copy operator: "
       << "poly1cp2 = " << poly1cp2.toString() << endl;
  cout << "Copy operator: "
       << "poly2cp2 = " << poly2cp2.toString() << endl;
  poly1cp2 *= op5;
  poly2cp2 *= op4;
  cout << "*= operator: (" << poly1.toString() << " *= " << op5.toString()
       << ") = " << poly1cp2.toString() << endl;
  cout << "*= operator: (" << poly2.toString() << " *= " << op4.toString()
       << ") = " << poly2cp2.toString() << endl;
  cout << "poly1cp2 = " << poly1cp2.toString() << endl;
  cout << "poly2cp2 = " << poly2cp2.toString() << endl;
  poly1cp2 *= poly2cp2;
  cout << "*= operator: (poly1cp2 *= poly2cp2) =  " << poly1cp2.toString() << endl;

  /*Normal Order tests.*/
  cout << "\nNormal Order Tests" << endl;
  cout << "poly1cp2 = " << poly1cp2.toString() << endl;
  poly1cp2.normalOrder();
  cout << "Normal Order poly1cp2 => " << poly1cp2.toString() << endl;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly11(mn1);
  mn1.herm();
  poly11 *= mn1;
  cout << "poly11 = " << poly11.toString() << endl;
  poly11.normalOrder();
  poly11.eraseNonNorm();
  cout << "poly11 => " << poly11.toString() << endl;
  poly1cp2 *= poly11;
  poly1cp2.eraseZeros();
  cout << "poly1cp2 *= poly11 = " << poly1cp2.toString() << endl;
  poly1cp2.normalOrder();
  cout << "poly1cp2 => " << poly1cp2.toString() << endl;

  /*Exit*/
  cout << "\nTests pass!" << endl;
  return EXIT_SUCCESS;
}

#endif  //QM_FERMI_POLYNOMIAL_TEST_CPP
