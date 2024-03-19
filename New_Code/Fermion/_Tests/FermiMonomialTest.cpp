#include "../fermiOperators.hpp"

using std::cout;
using std::endl;

int main(void) {
  Fermi1DLadderOp op1(1, true);
  Fermi1DLadderOp op2(2, false);
  Fermi1DLadderOp op3(3, false);
  Fermi1DLadderOp op4(4, true);
  Fermi1DLadderOp op5(5, true);
  /*Constructors tests*/
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
  /*Overloaded operators tests*/
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
  cout << "Tests pass!" << endl;
  mn1cp.herm();
  mn2cp.herm();
  cout << "herm operator: "
       << "mn1cp.herm() = " << mn1cp.toString() << endl;
  cout << "herm operator: "
       << "mn2cp.herm() = " << mn2cp.toString() << endl;
  return EXIT_SUCCESS;
}
