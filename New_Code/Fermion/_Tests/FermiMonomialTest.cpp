#include "../fermiOperators.hpp"

using std::cout;
using std::endl;

int main(void) {
  Fermi1DLadderOp op1(1, true);
  Fermi1DLadderOp op2(2, false);
  Fermi1DLadderOp op3(3, false);
  Fermi1DLadderOp op4(4, true);
  Fermi1DLadderOp op5(5, true);
  /*Constructors tests.*/
  FermiMonomial<Fermi1DLadderOp> mn0;
  FermiMonomial<Fermi1DLadderOp> mn1(op1);
  FermiMonomial<Fermi1DLadderOp> mn2(op2);
  cout << "Constructors tests" << endl;
  cout << "Default constructor: "
       << "mn0 = " << mn0.toString() << endl;
  cout << "Explicit constructor: "
       << "mn1 = " << mn1.toString() << endl;
  cout << "Explicit constructor: "
       << "mn2 = " << mn2.toString() << endl;
  FermiMonomial<Fermi1DLadderOp> mn1cp(mn1);
  FermiMonomial<Fermi1DLadderOp> mn2cp(mn2);
  cout << "Copy constructor: "
       << "mn1cp = " << mn1cp.toString() << endl;
  cout << "Copy constructor: "
       << "mn2cp = " << mn2cp.toString() << endl;

  /*Overloaded operators tests.*/
  FermiMonomial<Fermi1DLadderOp> mn1cp2 = mn1;
  FermiMonomial<Fermi1DLadderOp> mn2cp2 = mn2;
  cout << "\nOverloaded operators tests" << endl;
  cout << "Copy operator: "
       << "mn1cp2 = " << mn1cp2.toString() << endl;
  cout << "Copy operator: "
       << "mn2cp2 = " << mn2cp2.toString() << endl;
  cout << "Equal operator: (" << mn1.toString() << " == " << mn1cp.toString()
       << ") = " << (mn1 == mn1cp) << endl;
  cout << "Equal operator: (" << mn1.toString() << " != " << mn1cp.toString()
       << ") = " << (mn1 != mn1cp) << endl;
  cout << "Equal operator: (" << mn1.toString() << " == " << mn2.toString()
       << ") = " << (mn1 == mn2) << endl;
  cout << "Equal operator: (" << mn1.toString() << " != " << mn2.toString()
       << ") = " << (mn1 != mn2) << endl;
  mn1cp *= op5;
  mn2cp *= mn1cp;
  cout << "*= operator: (" << mn1.toString() << " *= " << op5.toString()
       << ") = " << mn1cp.toString() << endl;
  cout << "*= operator: (" << mn2.toString() << " *= " << mn1cp.toString()
       << ") = " << mn2cp.toString() << endl;
  mn1cp.herm();
  mn2cp.herm();
  cout << "herm operator: "
       << "mn1cp.herm() = " << mn1cp.toString() << endl;
  cout << "herm operator: "
       << "mn2cp.herm() = " << mn2cp.toString() << endl;

  /*Normalization tests.*/
  mn1cp.herm();
  mn2cp.herm();
  cout << "\nNormalization tests" << endl;
  cout << "Find wrong order: (" << mn1cp.toString()
       << ".FindWrongorder() = " << mn1cp.findWrongOrder()
       << "), isNorm = " << mn1cp.isNorm() << endl;
  cout << "Find wrong order: (" << mn2cp.toString()
       << ".FindWrongorder() = " << mn2cp.findWrongOrder()
       << "), isNorm = " << mn2cp.isNorm() << endl;
  mn1cp.herm();
  mn2cp.herm();
  cout << "Find wrong order: (" << mn1cp.toString()
       << ".FindWrongorder() = " << mn1cp.findWrongOrder()
       << "), isNorm = " << mn1cp.isNorm() << endl;
  cout << "Find wrong order: (" << mn2cp.toString()
       << ".FindWrongorder() = " << mn2cp.findWrongOrder()
       << "), isNorm = " << mn2cp.isNorm() << endl;
  cout << "Slice: "
       << "start to 1 " << mn2cp.sliceExprS(1).toString() << endl;
  cout << "Slice: "
       << "start to 2 " << mn2cp.sliceExprS(2).toString() << endl;
  cout << "Slice: "
       << "start to 3 " << mn2cp.sliceExprS(3).toString() << endl;
  cout << "Slice: "
       << "0 to end " << mn2cp.sliceExprE(0).toString() << endl;
  cout << "Slice: "
       << "1 to end " << mn2cp.sliceExprE(1).toString() << endl;
  cout << "Slice: "
       << "2 to end " << mn2cp.sliceExprE(2).toString() << endl;
  FermiMonomial<Fermi1DLadderOp> test1(mn2cp.sliceExprS(2));
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly1(complex<double>(1, 0), test1);
  cout << "poly1 = " << poly1.toString() << endl;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly2 = poly1;
  cout << "poly2 = " << poly2.toString() << endl;

  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly01 =
      NormOnce(complex<double>(1, 0), mn1cp);
  cout << "Norm once: " << mn1cp.toString() << " => "
       << "poly1 = " << poly01.toString() << endl;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly02 =
      NormOnce(complex<double>(1, 0), mn2cp);
  poly02.eraseZeros();
  cout << "Norm once: " << mn2cp.toString() << " => "
       << "poly1 = " << poly02.toString() << endl;
  cout << "\nTests pass!" << endl;
  return EXIT_SUCCESS;
}
