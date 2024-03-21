#include "../fermiStates.hpp"

using std::cout;
using std::endl;

int main(void) {
  /*Constructors tests*/
  FermiFockState state0;
  vector<bool> vec1;
  vec1.push_back(0);
  vec1.push_back(1);
  vec1.push_back(1);
  vector<bool> vec2(vec1);
  vec2[1] = 0;
  FermiFockState state1(vec1);
  FermiFockState state2(vec2);
  FermiFockState state1cp(state1);
  FermiFockState state2cp(state2);
  cout << "Constructors tests" << endl;
  cout << "Default constructor: "
       << "state0 = " << state0.toString() << endl;
  cout << "Explicit constructor: "
       << "state1 = " << state1.toString() << endl;
  cout << "Explicit constructor: "
       << "state2 = " << state2.toString() << endl;
  cout << "Copy constructor: "
       << "state1cp = " << state1cp.toString() << endl;
  cout << "Copy constructor: "
       << "state2cp = " << state2cp.toString() << endl;
  /*Overloaded operators tests*/
  FermiFockState state1cp2 = state1;
  FermiFockState state2cp2 = state2;
  cout << "\nOverloaded operators test" << endl;
  cout << "Copy operator: "
       << "state1cp2 = " << state1cp2.toString() << endl;
  cout << "Copy operator: "
       << "state2cp2 = " << state2cp2.toString() << endl;
  cout << "Equal operator: (" << state1.toString() << " == " << state1cp.toString()
       << ") = " << (state1 == state1cp) << endl;
  cout << "Equal operator: (" << state1.toString() << " != " << state1cp.toString()
       << ") = " << (state1 != state1cp) << endl;
  cout << "Equal operator: (" << state1.toString() << " == " << state2cp.toString()
       << ") = " << (state1 == state2cp) << endl;
  cout << "Equal operator: (" << state1.toString() << " != " << state2cp.toString()
       << ") = " << (state1 != state2cp) << endl;
  cout << "Tests pass!" << endl;
  return EXIT_SUCCESS;
}
