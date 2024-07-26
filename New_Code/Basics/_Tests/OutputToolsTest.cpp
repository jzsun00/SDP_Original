/*
  Jiazheng Sun
  Updated: Jul 26, 2024
  
  Test output tools, print complex, double and int numbers, vectors and matrices.
*/

#ifndef ORI_SDP_GS_SETTINGS_OUTPUTTOOLSTEST_CPP
#define ORI_SDP_GS_SETTINGS_OUTPUTTOOLSTEST_CPP

#include "../settings.hpp"

using std::complex;
using std::cout;
using std::endl;
using std::vector;

int main(void) {
  cout << "Test Output Tools" << endl;
  cout << "COMPLEX_PRECISION = " << COMPLEX_PRECISION << endl << endl;

  //Print single complex<double> numbers
  complex<double> num1(1, 2);
  complex<double> num2(2, -3);
  complex<double> num3(-3, 4);
  complex<double> num4(-4, -5);
  cout << "Print single complex<double> number" << endl;
  cout << "num1(++) = " << complex_toString(num1) << endl;
  cout << "num2(+-) = " << complex_toString(num2) << endl;
  cout << "num3(-+) = " << complex_toString(num3) << endl;
  cout << "num4(--) = " << complex_toString(num4) << endl << endl;

  //Print int vectors
  cout << "Print int vector" << endl;
  vector<int> intVec1{0, 1, 2, 3, 4, 5};
  vector<int> intVec2{0, 1, -2, 3, -4, 5};
  cout << "intVec1(+) =\n" << intVector_toString(intVec1) << endl;
  cout << "intVec2(+-) =\n" << intVector_toString(intVec2) << endl << endl;

  //Print double vectors
  cout << "Print double vector" << endl;
  vector<double> doubleVec1{0, 1.2, 2.3, 3.4, 4.5, 5.5};
  vector<double> doubleVec2{0, -1.2, 2.3, -3.4, 4.5, -5.6};
  cout << "doubleVec1(+) =\n" << doubleVector_toString(doubleVec1) << endl;
  cout << "doubleVec2(+-) =\n" << doubleVector_toString(doubleVec2) << endl << endl;

  //Print complex<double> vectors
  cout << "Print complex<double> vector";
  vector<complex<double> > vec1;
  vec1.push_back(num1);
  vec1.push_back(num2);
  vector<complex<double> > vec2;
  vec2.push_back(num3);
  vec2.push_back(num4);
  cout << "complexVec1(num1, num2) =\n" << complexVector_toString(vec1) << endl;
  cout << "complexVec2(num3, num4) =\n" << complexVector_toString(vec2) << endl << endl;

  //Print complex<double> matrices
  vector<vector<complex<double> > > mat1;
  mat1.push_back(vec1);
  vector<vector<complex<double> > > mat2;
  mat2.push_back(vec1);
  mat2.push_back(vec2);
  vector<vector<complex<double> > > mat3;
  mat3.push_back(vec1);
  mat3.push_back(vec2);
  mat3.push_back(vec1);
  cout << "mat1(vec1) =\n" << complexMatrix_toString(mat1) << endl << endl;
  cout << "mat2(vec1, vec2) =\n" << complexMatrix_toString(mat2) << endl << endl;
  cout << "mat3(vec1, vec2, vec1) =\n" << complexMatrix_toString(mat3) << endl << endl;

  //Exit
  cout << "Tests pass!" << endl;
  return EXIT_SUCCESS;
}

#endif  //ORI_SDP_GS_SETTINGS_OUTPUTTOOLSTEST_CPP
