/*
  Jiazheng Sun
  Updated: Jun 26, 2024

  Class Implementations:
  SparseHamiltonian<PolyType, BaseStateType>
  FullHamiltonian.<PolyType, BaseStateType>
*/

#ifndef QM_HAMILTONIANS_TEM_HPP
#define QM_HAMILTONIANS_TEM_HPP

#include "./hamiltonians.hpp"

//------------------------------------SparseHamiltonian<PolyType, BaseStateType>---------

template<typename PolyType, typename BaseStateType>
std::string SparseHamiltonian<PolyType, BaseStateType>::toString() {
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

template<typename polyType, typename BaseStateType>
void SparseHamiltonian<polyType, BaseStateType>::createMatrix(
    Basis<BaseStateType> & basis) {
  dim = basis.getSize();
  pcol.push_back(0);
  for (long unsigned j = 0; j < dim; j++) {
    State<BaseStateType> mid = poly * basis[j];
    for (long unsigned i = 0; i < dim; i++) {
      std::complex<double> elementij = innerProduct(basis[i], mid);
      if (std::abs(elementij) <= ERROR) {
        continue;
      }
      nnz++;
      nzVal.push_back(elementij);
      irow.push_back(i);
    }
    pcol.push_back(nnz);
  }
  //pcol.push_back(nnz);
}

//----------------------------------SparseRealHamiltonian<PolyType, BaseStateType>-------

template<typename polyType, typename BaseStateType>
void SparseRealHamiltonian<polyType, BaseStateType>::createMatrix(
    Basis<BaseStateType> & basis) {
  dim = basis.getSize();
  pcol.push_back(0);
  for (long unsigned j = 0; j < dim; j++) {
    State<BaseStateType> mid = poly * basis[j];
    for (long unsigned i = 0; i < dim; i++) {
      std::complex<double> elementij = innerProduct(basis[i], mid);
      if (std::abs(elementij) <= ERROR) {
        continue;
      }
      nnz++;
      nzVal.push_back(elementij.real());
      irow.push_back(i);
    }
    pcol.push_back(nnz);
  }
}

//-------------------------------------FullHamiltonian<PolyType, BaseStateType>----------

template<typename PolyType, typename BaseStateType>
std::string FullHamiltonian<PolyType, BaseStateType>::toString() {
  std::string ans;
  ans += "The Polynomial:\n";
  ans += poly.toString();
  ans += "\nSparsity: ";
  ans += std::to_string((double)nnz / (dim * dim));
  ans += "\n";
  ans += "The Full Matrix:\n";
  for (long unsigned i = 0; i < dim; i++) {
    for (long unsigned j = 0; j < dim; j++) {
      ans += "  " + std::to_string(matrix[i][j].real()) + "+i" +
             std::to_string(matrix[i][j].imag());
    }
    ans += "\n";
  }
  return ans;
}

template<typename PolyType, typename BaseStateType>
void FullHamiltonian<PolyType, BaseStateType>::createMatrix(
    Basis<BaseStateType> const & basis) {
  for (long unsigned i = 0; i < dim; i++) {
    for (long unsigned j = 0; j < dim; j++) {
      //std::cout << "i = " << i << ", j = " << j << std::endl;
      matrix[i][j] = innerProduct(basis[i], poly * basis[j]);
      if (std::abs(matrix[i][j]) > ERROR) {
        nnz++;
      }
    }
  }
}

#endif  //QM_HAMILTONIANS_TEM_HPP
