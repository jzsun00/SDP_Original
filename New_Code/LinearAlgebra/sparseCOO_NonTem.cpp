/*
  Jiazheng Sun
  Updated: Jul 22, 2024

  Class Implementations:
  ComplexCOOMatrix
*/

#ifndef LA_SPARSE_MATRIX_COO_NON_TEM_CPP
#define LA_SPARSE_MATRIX_COO_NON_TEM_CPP

#include "./printTools.hpp"
#include "./sparseCOO.hpp"

//-----------------------------------------------------------ComplexCOOMatrix------------

std::string ComplexCOOMatrix::element_toString(std::complex<double> element) const {
  return LA::complex_toString(element);
}

#endif  //LA_SPARSE_MATRIX_COO_NON_TEM_CPP
