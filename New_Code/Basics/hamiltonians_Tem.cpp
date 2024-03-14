#ifndef ORI_SDP_GS_HAMILTONIANS_CPP
#define ORI_SDP_GS_HAMILTONIANS_CPP

#include "./hamiltonians.hpp"

//-------------------------------------------------------------SparseHamiltonian---------
/*
template<typename polyType, typename basisType>
void SparseHamiltonian<polyType, basisType>::createMatrix(basisType & basis) {
  double error = std::pow(10, -12);
  for (long unsigned i = 0; i < dim; i++) {
    for (long unsigned j = 0; j < dim; j++) {
      complex<double> element = innerProduct(basis[i], poly * basis[j]);
      if (std::abs(element) > error) {
        nnz++;
      }
    }
  }
}
*/

//-------------------------------------------------------------FullHamiltonian-----------

template<typename polyType, typename basisType>
std::string FullHamiltonian<polyType, basisType>::toString() {
  std::string ans = "";
  ans += "The Full Matrix:\n";
  for (long unsigned i = 0; i < dim; i++) {
    for (long unsigned j = 0; j < dim; j++) {
      ans += "  " + std::to_string(matrix[i][j].real()) + "+" +
             std::to_string(matrix[i][j].imag());
    }
    ans += "\n";
  }
  return ans;
}

template<typename polyType, typename basisType>
void FullHamiltonian<polyType, basisType>::createMatrix(basisType const & basis) {
  for (long unsigned i = 0; i < dim; i++) {
    for (long unsigned j = 0; j < dim; j++) {
      //std::cout << "i = " << i << ", j = " << j << std::endl;
      matrix[i][j] = innerProduct(basis[i], poly * basis[j]);
    }
  }
}

#endif
