/*
  Jiazheng Sun
  Updated: Aug 6, 2024
  
  Test implementations of FermiMonomial<OpType>.
*/

#ifndef QM_FERMI_MONOMIAL_TEST_CPP
#define QM_FERMI_MONOMIAL_TEST_CPP

#include <exception>

#include "../fermiOperators_Tem.hpp"

using std::complex;
using std::cout;
using std::endl;

int main(void) {
  cout << "Fermi Monomial Test\n" << endl;

  /*Constructors Tests.*/
  Fermi1DLadderOp op1(1, true);
  Fermi1DLadderOp op2(2, false);
  Fermi1DLadderOp op3(3, false);
  Fermi1DLadderOp op4(4, true);
  Fermi1DLadderOp op5(5, true);
  FermiMonomial<Fermi1DLadderOp> mn0;
  FermiMonomial<Fermi1DLadderOp> mn1(op1);
  FermiMonomial<Fermi1DLadderOp> mn2(op2);
  cout << "Constructors Tests" << endl;
  cout << "Default constructor: "
       << "mn0 = " << mn0.toString() << endl;
  cout << "Explicit constructor(1, true): "
       << "mn1 = " << mn1.toString() << endl;
  cout << "Explicit constructor(2, false): "
       << "mn2 = " << mn2.toString() << endl;
  FermiMonomial<Fermi1DLadderOp> mn1cp(mn1);
  FermiMonomial<Fermi1DLadderOp> mn2cp(mn2);
  cout << "Copy constructor: "
       << "mn1cp(mn1) = " << mn1cp.toString() << endl;
  cout << "Copy constructor: "
       << "mn2cp(mn2) = " << mn2cp.toString() << endl;

  /*Overloaded Operators Tests.*/
  FermiMonomial<Fermi1DLadderOp> mn1cp2 = mn1;
  FermiMonomial<Fermi1DLadderOp> mn2cp2 = mn2;
  cout << "\nOverloaded Operators Tests" << endl;
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
  cout << "mn1cp = " << mn1cp.toString() << endl;
  cout << "mn2cp = " << mn2cp.toString() << endl;
  mn1cp.herm();
  mn2cp.herm();
  cout << "herm operator: "
       << "mn1cp.herm() = " << mn1cp.toString() << endl;
  cout << "herm operator: "
       << "mn2cp.herm() = " << mn2cp.toString() << endl;

  /*Normal Order Tests.*/
  cout << "\nNormal Order Tests" << endl;
  mn1cp.herm();
  mn2cp.herm();
  mn1cp.reverse();
  mn2cp.reverse();
  cout << "reverse operator: "
       << "mn1cp.reverse() = " << mn1cp.toString() << endl;
  cout << "reverse operator: "
       << "mn2cp.reverse() = " << mn2cp.toString() << endl;
  mn1cp.reverse();
  mn2cp.reverse();
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
  cout << "mn2cp = " << mn2cp.toString() << endl;
  cout << "Slice: "
       << "start to 1 " << mn2cp.sliceExprStart(1).toString() << endl;
  cout << "Slice: "
       << "start to 2 " << mn2cp.sliceExprStart(2).toString() << endl;
  cout << "Slice: "
       << "start to 3 " << mn2cp.sliceExprStart(3).toString() << endl;
  cout << "Slice: "
       << "0 to end " << mn2cp.sliceExprEnd(0).toString() << endl;
  cout << "Slice: "
       << "1 to end " << mn2cp.sliceExprEnd(1).toString() << endl;
  cout << "Slice: "
       << "2 to end " << mn2cp.sliceExprEnd(2).toString() << endl;
  try {
    FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly01 =
        NormOnce(complex<double>(1.0, 0), mn1cp);
    cout << "Norm once: " << mn1cp.toString() << " => " << poly01.toString() << endl;
    FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly02 =
        NormOnce(complex<double>(1.0, 0), mn2cp);
    cout << "Norm once: " << mn2cp.toString() << " => " << poly02.toString() << endl;
    mn1cp.reverse();
    FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly03 =
        NormOnce(complex<double>(1.0, 0), mn1cp);
    cout << "Norm once: " << mn1cp.toString() << " => " << poly03.toString() << endl;
  }
  catch (std::exception & e) {
    cout << e.what() << endl;
  }

  /*Infinite System Tests.*/
  cout << "\nInfinite System Tests" << endl;

  /*Exit*/
  cout << "\nTests pass!" << endl;
  return EXIT_SUCCESS;
}

#endif  //QM_FERMI_MONOMIAL_TEST_CPP
