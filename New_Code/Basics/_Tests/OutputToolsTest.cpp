/*
  Jiazheng Sun
  Updated: Jun 10, 2024
  
  Test output tools, print complex numbers, vectors and matrices.
*/

#ifndef ORI_SDP_GS_OUTPUTTOOLSTEST_CPP
#define ORI_SDP_GS_OUTPUTTOOLSTEST_CPP

#include "../settings.hpp"

using std::cout;
using std::endl;

int main(void) {
  cout << "Test Output Tools" << endl;
  cout << "COMPLEX_PRECISION = " << COMPLEX_PRECISION << endl << endl;
  //Print single complex<double> numbers
  complex<double> num1(1, 2);
  complex<double> num2(2, -3);
  complex<double> num3(-3, 4);
  complex<double> num4(-4, -5);
  cout << "num1(++) = " << complex_toString(num1) << endl;
  cout << "num2(+-) = " << complex_toString(num2) << endl;
  cout << "num3(-+) = " << complex_toString(num3) << endl;
  cout << "num4(--) = " << complex_toString(num4) << endl << endl;
  //Print complex<double> vectors
  vector<complex<double> > vec1;
  vec1.push_back(num1);
  vec1.push_back(num2);
  vector<complex<double> > vec2;
  vec2.push_back(num3);
  vec2.push_back(num4);
  cout << "vec1(num1, num2) = " << endl << complexVector_toString(vec1) << endl;
  cout << "vec2(num3, num4) = " << endl << complexVector_toString(vec2) << endl << endl;
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
  cout << "mat1(vec1) = " << endl << complexMatrix_toString(mat1) << endl << endl;
  cout << "mat2(vec1, vec2) = " << endl << complexMatrix_toString(mat2) << endl << endl;
  cout << "mat3(vec1, vec2, vec1) = " << endl
       << complexMatrix_toString(mat3) << endl
       << endl;
  //Exit
  cout << "Tests pass!" << endl;
  return EXIT_SUCCESS;
}

#endif  //ORI_SDP_GS_OUTPUTTOOLSTEST_CPP
