#ifndef ORI_SDP_GS_HAMILTONIANS_CPP
#define ORI_SDP_GS_HAMILTONIANS_CPP

#include "./hamiltonians.hpp"

template<typename polyType, typename basisType>
void Hamiltonian<polyType, basisType>::createMatrix(basisType & basis) {
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

#endif
