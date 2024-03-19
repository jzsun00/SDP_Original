#include <cassert>

#include "../fermiOperators.hpp"

using std::cout;
using std::endl;

int main(void) {
  /*Constructors tests.*/
  Fermi1DLadderOp lad;
  Fermi1DLadderOp lad1(1, true);
  Fermi1DLadderOp lad2(2, false);
  Fermi1DLadderOp lad1cp(lad1);
  Fermi1DLadderOp lad2cp(lad2);
  cout << "Constructors test" << endl;
  cout << "Default constructor: "
       << "lad = " << lad.toString() << endl;
  cout << "Explicit constructor: "
       << "lad1 = " << lad1.toString() << endl;
  cout << "Explicit constructor: "
       << "lad2 = " << lad2.toString() << endl;
  cout << "Copy constructor: "
       << "lad1cp = " << lad1cp.toString() << endl;
  cout << "Copy constructor: "
       << "lad2cp = " << lad2cp.toString() << endl;
  /*Overloaded operators test.*/
  Fermi1DLadderOp lad1cp2 = lad1;
  Fermi1DLadderOp lad2cp2 = lad2;
  Fermi1DLadderOp lad3(3, false);
  Fermi1DLadderOp lad4(4, true);
  cout << "\nOverloaded operators test" << endl;
  cout << "Copy operator: "
       << "lad1cp2 = " << lad1cp2.toString() << endl;
  cout << "Copy operator: "
       << "lad2cp2 = " << lad2cp2.toString() << endl;
  lad1cp.herm();
  lad2cp.herm();
  cout << "Equal operator: (" << lad1.toString() << " == " << lad1cp.toString()
       << ") = " << (lad1 == lad1cp) << endl;
  cout << "Equal operator: (" << lad1.toString() << " != " << lad1cp.toString()
       << ") = " << (lad1 != lad1cp) << endl;
  cout << "Equal operator: (" << lad1.toString() << " == " << lad2cp.toString()
       << ") = " << (lad1 == lad2cp) << endl;
  cout << "Equal operator: (" << lad1.toString() << " != " << lad2cp.toString()
       << ") = " << (lad1 != lad2cp) << endl;
  cout << "Equal operator: (" << lad1.toString() << " == " << lad1cp2.toString()
       << ") = " << (lad1 == lad1cp2) << endl;
  cout << "Equal operator: (" << lad1.toString() << " != " << lad1cp2.toString()
       << ") = " << (lad1 != lad1cp2) << endl;
  cout << "Less than operator: (" << lad1.toString() << " < " << lad1cp.toString()
       << ") = " << (lad1 < lad1cp) << endl;
  cout << "Less than operator: (" << lad1.toString() << " > " << lad1cp.toString()
       << ") = " << (lad1 > lad1cp) << endl;
  cout << "Less than operator: (" << lad1.toString() << " < " << lad2cp.toString()
       << ") = " << (lad1 < lad2cp) << endl;
  cout << "Less than operator: (" << lad1.toString() << " > " << lad2cp.toString()
       << ") = " << (lad1 > lad2cp) << endl;
  cout << "Less than operator: (" << lad2cp.toString() << " < " << lad3.toString()
       << ") = " << (lad2cp < lad3) << endl;
  cout << "Less than operator: (" << lad2cp.toString() << " > " << lad3.toString()
       << ") = " << (lad2cp > lad3) << endl;
  cout << "Tests pass!" << endl;
  return EXIT_SUCCESS;
}
