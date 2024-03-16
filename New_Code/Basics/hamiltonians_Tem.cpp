#ifndef ORI_SDP_GS_HAMILTONIANS_CPP
#define ORI_SDP_GS_HAMILTONIANS_CPP

#include "./hamiltonians.hpp"

//-------------------------------------------------------------SparseHamiltonian---------

template<typename PolyType, typename BasisType>
std::string SparseHamiltonian<PolyType, BasisType>::toString() {
  std::string ans = "";
  ans += "Number of non-zero elements: " + std::to_string(nnz) + "\n";
  ans += "Non-zero values:\n";
  for (size_t i = 0; i < nzVal.size(); i++) {
    ans += "  " + std::to_string(nzVal[i].real()) + " + " +
           std::to_string(nzVal[i].imag()) + "i";
  }
  ans += "\nRow indices of non-zero elements:\n";
  for (size_t i = 0; i < irow.size(); i++) {
    ans += "  " + std::to_string(irow[i]);
  }
  ans += "\nBeginning of each coloumn:\n";
  for (size_t i = 0; i < pcol.size(); i++) {
    ans += "  " + std::to_string(pcol[i]);
  }
  return ans;
}

template<typename polyType, typename basisType>
void SparseHamiltonian<polyType, basisType>::createMatrix(basisType & basis) {
  pcol.push_back(0);
  for (long unsigned i = 0; i < dim; i++) {
    for (long unsigned j = 0; j < dim; j++) {
      complex<double> elementij = innerProduct(basis[i], poly * basis[j]);
      if (std::abs(elementij) <= ERROR) {
        continue;
      }
      nnz++;
      nzVal.push_back(elementij);
      irow.push_back(j);
    }
    pcol.push_back(nnz);
  }
  //pcol.push_back(nnz);
}

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
