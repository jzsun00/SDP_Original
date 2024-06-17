/*
  Jiazheng Sun
  Updated: Jun 17, 2024

  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
#define ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP

#include <cstddef>
#include <map>
#include <vector>

#include "./hamiltonians_XXZ.hpp"

//------------------------------------------------------------XXZSparseHamiltonian-------

void XXZSparseHamiltonian::createMatrix(SpinHalfBasis1D & basis) {
  dim = basis.getSize();
  pcol.push_back(0);
  for (long unsigned j = 0; j < dim; j++) {
    SpinHalfState1D mid = poly * basis[j];
    //mid.eraseZeros();

    std::vector<int> indices;
    std::map<int, complex<double> > elements;
    for (size_t k = 0; k < mid.getSize(); k++) {
      int index = basis.findBaseState(mid[k].second);
      complex<double> value = mid[k].first;
      indices.push_back(index);
      elements[index] = value;
    }
    std::sort(indices.begin(), indices.end());
    for (size_t k = 0; k < mid.getSize(); k++) {
      nnz++;
      nzVal.push_back(elements[indices[k]]);
      irow.push_back(indices[k]);
    }

    /*
    for (long unsigned i = 0; i < dim; i++) {
      complex<double> elementij = innerProduct(basis[i], mid);
      if (std::abs(elementij) <= ERROR) {
        continue;
      }
      nnz++;
      nzVal.push_back(elementij);
      irow.push_back(i);
    }
    */
    pcol.push_back(nnz);
  }
  //pcol.push_back(nnz);
}

//----------------------------------------------------------------Other Functions--------

SpinHalfPolynomial1D makePoly(size_t sites, double Jz) {
  SpinHalfPolynomial1D ans;
  for (size_t i = 0; i < sites - 1; i++) {
    SpinHalfOp1D Sz(i);
    SpinHalfOp1D SzN(i + 1);
    SpinHalfOp1D Su(i, true);
    SpinHalfOp1D Sd(i, false);
    SpinHalfOp1D SuN(i + 1, true);
    SpinHalfOp1D SdN(i + 1, false);
    SpinHalfMonomial1D MNud(Su);
    MNud *= SdN;
    SpinHalfMonomial1D MNdu(Sd);
    MNdu *= SuN;
    SpinHalfMonomial1D MNz(Sz);
    MNz *= SzN;
    if (i == 0 || i == sites - 2) {
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.25, 0), MNud);
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.25, 0), MNdu);
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(Jz / 2, 0), MNz);
    }
    else {
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.5, 0), MNud);
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.5, 0), MNdu);
      ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(Jz, 0), MNz);
    }
  }
  return ans;
}

#endif  //ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
