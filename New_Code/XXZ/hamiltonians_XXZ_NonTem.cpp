/*
  Jiazheng Sun
  Updated: Aug 16, 2024
  
  Class Implementations:
  XXZSparseHamiltonian
  XXZSparseRealHamiltonian
  XXZFullHamiltonian
  
  Function Implementations:
  SpinHalfPolynomial1D makeSpinPoly(size_t sites, double Jz);
  SpinHalfPolynomial1D makeSpinPolyPBC(size_t sites, double Jz);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > makeFermiPoly(int start, int end, double Jz);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > makeFermiFinitePoly(int start, int end, double Jz);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > makeFermiPolyPBC(int start, int end, double Jz);
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > makeHardCorePoly(size_t sites, double Jz);
  SpinHalfState1D makeMidState(size_t sites, double Jz, const SpinHalfBaseState1D & rhs);
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
#define ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP

#include "./hamiltonians_XXZ.hpp"
#include "omp.h"

using std::complex;
using std::pair;
using std::vector;

//------------------------------------------------------------XXZSparseHamiltonian-------

void XXZSparseHamiltonian::createMatrix(SpinHalfBasis1D & basis) {
  dim = basis.getSize();
  size_t batchSize = 200000;
  size_t batchNum = dim / batchSize + 1;
  pcol.push_back(0);

  for (size_t batchIdx = 0; batchIdx < batchNum; batchIdx++) {
    size_t currentSize = batchSize;
    if (batchIdx == batchNum - 1) {
      currentSize = dim % batchSize;
    }
    vector<SpinHalfState1D> midStates(currentSize);

#pragma omp parallel
    {
#pragma omp for
      for (long unsigned j = 0; j < currentSize; j++) {
        midStates[j] = XXZ1D::makeMidState(sites, Jz, basis[j + batchSize * batchIdx]);
      }
    }

    for (long unsigned j = 0; j < currentSize; j++) {
      std::vector<int> indices;
      std::map<int, complex<double> > elements;
      for (size_t k = 0; k < midStates[j].getSize(); k++) {
        //int index = basis.findBaseState(mid[k].second);
        size_t index = basis.lookUpBaseState(midStates[j][k].second);
        complex<double> value = midStates[j][k].first;
        indices.push_back(index);
        elements[index] = value;
      }
      std::sort(indices.begin(), indices.end());
      for (size_t k = 0; k < midStates[j].getSize(); k++) {
        nnz++;
        nzVal.push_back(elements[indices[k]]);
        irow.push_back(indices[k]);
      }
      pcol.push_back(nnz);
    }
  }
  //pcol.push_back(nnz);
}

//------------------------------------------------------XXZSparseRealHamiltonian---------

void XXZSparseRealHamiltonian::createMatrix(SpinHalfBasis1D & basis) {
  dim = basis.getSize();
  size_t batchSize = 500000;
  size_t batchNum = dim / batchSize + 1;
  pcol.push_back(0);

  for (size_t batchIdx = 0; batchIdx < batchNum; batchIdx++) {
    size_t currentSize = batchSize;
    if (batchIdx == batchNum - 1) {
      currentSize = dim % batchSize;
    }
    vector<SpinHalfState1D> midStates(currentSize);

//omp_set_num_threads(8);
#pragma omp parallel
    {
#pragma omp for
      for (long unsigned j = 0; j < currentSize; j++) {
        //midStates[j] = makeMidState(sites, Jz, basis[j + batchSize * batchIdx]);
        midStates[j] = poly * basis[j + batchSize * batchIdx];
      }
    }

    vector<vector<int> > indices(currentSize);
    vector<std::map<int, double> > elements(currentSize);

#pragma omp parallel
    {
#pragma omp for
      for (long unsigned j = 0; j < currentSize; j++) {
        //std::vector<int> indices;
        //std::map<int, double> elements;
        size_t midStateSize = midStates[j].getSize();
        for (size_t k = 0; k < midStateSize; k++) {
          //int index = basis.findBaseState(mid[k].second);
          size_t index = basis.lookUpBaseState(midStates[j][k].second);
          if (j + batchSize * batchIdx <= index) {
            complex<double> value = midStates[j][k].first;
            indices[j].push_back(index);
            elements[j][index] = value.real();
          }
        }
      }
    }

    for (long unsigned j = 0; j < currentSize; j++) {
      std::sort(indices[j].begin(), indices[j].end());
      size_t indicesSize = indices[j].size();
      for (size_t k = 0; k < indicesSize; k++) {
        nnz++;
        nzVal.push_back(elements[j][indices[j][k]]);
        irow.push_back(indices[j][k]);
      }
      pcol.push_back(nnz);
    }
  }
}

void XXZSparseRealHamiltonian::createFullBasisMatrix(const SpinHalfBasis1D & basis) {
  dim = basis.getSize();
  pcol.push_back(0);
  vector<SpinHalfState1D> midStates(dim);

#pragma omp parallel for
  for (size_t i = 0; i < dim; i++) {  //Parallel compute all (poly * baseState)
    midStates[i] = poly * basis[i];
  }

  for (size_t i = 0; i < dim; i++) {
    size_t colIdx = i;
    size_t midStateSize = midStates[i].getSize();
    vector<size_t> indices;
    std::unordered_map<size_t, double> elements;
    for (size_t j = 0; j < midStateSize; j++) {  //Fill in indices and elements
      size_t rowIdx = midStates[i][j].second.toDecimal();
      if (colIdx <= rowIdx) {
        indices.push_back(rowIdx);
        elements[rowIdx] = midStates[i][j].first.real();
      }
    }
    std::sort(indices.begin(), indices.end());
    size_t indicesSize = indices.size();
    for (size_t k = 0; k < indicesSize; k++) {  //Print into the matrix
      nnz++;
      nzVal.push_back(elements[indices[k]]);
      irow.push_back(indices[k]);
    }
    pcol.push_back(nnz);
  }
}

bool isReverseEqual(const SpinHalfBaseState1D & state) {
  size_t n = state.getSize();
  for (size_t i = 0; i < n / 2; ++i) {
    if (state[i] != state[n - 1 - i]) {
      return false;
    }
  }
  return true;
}

bool isLessOrEqualToReverse(const SpinHalfBaseState1D & state) {
  size_t n = state.getSize();
  for (size_t i = 0; i < n; ++i) {
    if (state[i] < state[n - 1 - i]) {
      return true;
    }
    else if (state[i] > state[n - 1 - i]) {
      return false;
    }
  }
  return true;  // They are equal
}

void XXZSparseRealHamiltonian::createRefSymMatrix(SpinHalfBasis1D & basis) {
  dim = basis.getSize();
  size_t batchSize = 500000;
  size_t batchNum = dim / batchSize + 1;
  pcol.push_back(0);

  for (size_t batchIdx = 0; batchIdx < batchNum; batchIdx++) {
    size_t currentSize = batchSize;
    if (batchIdx == batchNum - 1) {
      currentSize = dim % batchSize;
    }
    vector<SpinHalfState1D> midStates(currentSize);

#pragma omp parallel
    {
#pragma omp for
      for (long unsigned j = 0; j < currentSize; j++) {
        midStates[j] = XXZ1D::makeMidState(sites, Jz, basis[j + batchSize * batchIdx]);
      }
    }

    vector<vector<int> > indices(currentSize);
    vector<std::map<int, double> > elements(currentSize);

    bool symState = false;

#pragma omp parallel
    {
#pragma omp for
      for (long unsigned j = 0; j < currentSize; j++) {
        size_t midStateSize = midStates[j].getSize();
        bool symState1 = isReverseEqual(basis[j + batchSize * batchIdx]);
        for (size_t k = 0; k < midStateSize; k++) {
          //int index = basis.findBaseState(mid[k].second);
          size_t index;
          if (!isLessOrEqualToReverse(midStates[j][k].second)) {
            continue;
          }
          index = basis.lookUpBaseState(midStates[j][k].second);
          //std::cout << "index = " << index << std::endl;
          symState = isReverseEqual(midStates[j][k].second);
          if (j + batchSize * batchIdx <= index) {
            complex<double> value = midStates[j][k].first;
            if (symState1) {
              if (!symState) {
                value *= std::pow(2, 0.5);
              }
            }
            else {
              if (symState) {
                value *= std::pow(2, 0.5);
              }
            }
            indices[j].push_back(index);
            elements[j][index] = value.real();
          }
        }
      }
    }

    for (long unsigned j = 0; j < currentSize; j++) {
      std::sort(indices[j].begin(), indices[j].end());
      size_t indicesSize = indices[j].size();
      for (size_t k = 0; k < indicesSize; k++) {
        nnz++;
        nzVal.push_back(elements[j][indices[j][k]]);
        irow.push_back(indices[j][k]);
      }
      pcol.push_back(nnz);
    }
  }
}

//----------------------------------------------------------------Other Functions--------

/*Open boundary condition, interaction at the boundary is divided by 2.*/
SpinHalfPolynomial1D XXZ1D::makeSpinPoly(size_t sites, double Jz) {
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

/*Periodic boundary condition.*/
SpinHalfPolynomial1D XXZ1D::makeSpinPolyPBC(size_t sites, double Jz) {
  SpinHalfPolynomial1D ans;
  for (size_t i = 0; i < sites; i++) {
    SpinHalfOp1D Sz(i);
    SpinHalfOp1D SzN((i + 1) % sites);
    SpinHalfOp1D Su(i, true);
    SpinHalfOp1D Sd(i, false);
    SpinHalfOp1D SuN((i + 1) % sites, true);
    SpinHalfOp1D SdN((i + 1) % sites, false);
    SpinHalfMonomial1D MNud(Su);
    MNud *= SdN;
    SpinHalfMonomial1D MNdu(Sd);
    MNdu *= SuN;
    SpinHalfMonomial1D MNz(Sz);
    MNz *= SzN;
    ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.5, 0), MNud);
    ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(0.5, 0), MNdu);
    ans += pair<complex<double>, SpinHalfMonomial1D>(complex<double>(Jz, 0), MNz);
  }
  return ans;
}

FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > XXZ1D::makeFermiPoly(int start,
                                                                      int end,
                                                                      double Jz) {
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > ans;
  for (int i = start; i < end; i++) {
    Fermi1DLadderOp Su(i, true);
    Fermi1DLadderOp Sd(i, false);
    Fermi1DLadderOp SuN(i + 1, true);
    Fermi1DLadderOp SdN(i + 1, false);
    FermiMonomial<Fermi1DLadderOp> MNud(Su);
    MNud *= SdN;
    FermiMonomial<Fermi1DLadderOp> MNdu(Sd);
    MNdu *= SuN;
    FermiMonomial<Fermi1DLadderOp> MNz(Su);
    MNz *= Sd;
    FermiMonomial<Fermi1DLadderOp> MNzN(SuN);
    MNzN *= SdN;
    FermiMonomial<Fermi1DLadderOp> MNzz(MNz);
    MNzz *= MNzN;
    Fermi1DLadderOp unitOp(true);
    FermiMonomial<Fermi1DLadderOp> unitMn(unitOp);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(complex<double>(0.5, 0),
                                                                  MNud);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
        complex<double>(-0.5, 0), MNdu);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(complex<double>(Jz, 0),
                                                                  MNzz);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
        complex<double>(-0.5 * Jz, 0), MNz);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
        complex<double>(-0.5 * Jz, 0), MNzN);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
        complex<double>(0.25 * Jz, 0), unitMn);
  }
  return ans;
}

/*Open boundary condition in [start, end],
  interaction at the boundary is divided by 2.*/
FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > XXZ1D::makeFermiFinitePoly(int start,
                                                                            int end,
                                                                            double Jz) {
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > ans;
  for (int i = start; i < end; i++) {
    Fermi1DLadderOp Su(i, true);
    Fermi1DLadderOp Sd(i, false);
    Fermi1DLadderOp SuN(i + 1, true);
    Fermi1DLadderOp SdN(i + 1, false);
    FermiMonomial<Fermi1DLadderOp> MNud(Su);
    MNud *= SdN;
    FermiMonomial<Fermi1DLadderOp> MNdu(Sd);
    MNdu *= SuN;
    FermiMonomial<Fermi1DLadderOp> MNz(Su);
    MNz *= Sd;
    FermiMonomial<Fermi1DLadderOp> MNzN(SuN);
    MNzN *= SdN;
    FermiMonomial<Fermi1DLadderOp> MNzz(MNz);
    MNzz *= MNzN;
    Fermi1DLadderOp unitOp(true);
    FermiMonomial<Fermi1DLadderOp> unitMn(unitOp);
    if (i == start || i == end - 1) {
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(0.25, 0), MNud);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(-0.25, 0), MNdu);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(0.5 * Jz, 0), MNzz);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(-0.25 * Jz, 0), MNz);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(-0.25 * Jz, 0), MNzN);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(0.125 * Jz, 0), unitMn);
    }
    else {
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(0.5, 0), MNud);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(-0.5, 0), MNdu);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(Jz, 0), MNzz);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(-0.5 * Jz, 0), MNz);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(-0.5 * Jz, 0), MNzN);
      ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
          complex<double>(0.25 * Jz, 0), unitMn);
    }
  }
  return ans;
}

FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > XXZ1D::makeFermiPolyPBC(int start,
                                                                         int end,
                                                                         double Jz) {
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > ans;
  for (int i = start; i <= end; i++) {
    Fermi1DLadderOp Su(i, true);
    Fermi1DLadderOp Sd(i, false);
    Fermi1DLadderOp SuN(i + 1, true);
    Fermi1DLadderOp SdN(i + 1, false);
    if (i == end) {
      SuN = Fermi1DLadderOp(start, true);
      SdN = Fermi1DLadderOp(start, false);
    }
    FermiMonomial<Fermi1DLadderOp> MNud(Su);
    MNud *= SdN;
    FermiMonomial<Fermi1DLadderOp> MNdu(Sd);
    MNdu *= SuN;
    FermiMonomial<Fermi1DLadderOp> MNz(Su);
    MNz *= Sd;
    FermiMonomial<Fermi1DLadderOp> MNzN(SuN);
    MNzN *= SdN;
    FermiMonomial<Fermi1DLadderOp> MNzz(MNz);
    MNzz *= MNzN;
    Fermi1DLadderOp unitOp(true);
    FermiMonomial<Fermi1DLadderOp> unitMn(unitOp);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(complex<double>(0.5, 0),
                                                                  MNud);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
        complex<double>(-0.5, 0), MNdu);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(complex<double>(Jz, 0),
                                                                  MNzz);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
        complex<double>(-0.5 * Jz, 0), MNz);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
        complex<double>(-0.5 * Jz, 0), MNzN);
    ans += pair<complex<double>, FermiMonomial<Fermi1DLadderOp> >(
        complex<double>(0.25 * Jz, 0), unitMn);
  }
  return ans;
}

HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > XXZ1D::makeHardCorePoly(
    size_t sites,
    double Jz) {
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > ans;
  for (int i = 0; i < (int)sites - 1; i++) {
    //SpinHalfOp Sz(i);
    //SpinHalfOp SzN(i + 1);
    HardCore1DLadderOp Su(i, true);
    HardCore1DLadderOp Sd(i, false);
    HardCore1DLadderOp SuN(i + 1, true);
    HardCore1DLadderOp SdN(i + 1, false);
    HardCoreMonomial<HardCore1DLadderOp> MNud(Su);
    MNud *= SdN;
    HardCoreMonomial<HardCore1DLadderOp> MNdu(Sd);
    MNdu *= SuN;
    HardCoreMonomial<HardCore1DLadderOp> MNz(Su);
    MNz *= Sd;
    HardCoreMonomial<HardCore1DLadderOp> MNzN(SuN);
    MNzN *= SdN;
    HardCoreMonomial<HardCore1DLadderOp> MNzz(MNz);
    MNzz *= MNzN;
    HardCore1DLadderOp unitOp(true);
    HardCoreMonomial<HardCore1DLadderOp> unitMn(unitOp);
    /*
    if (i == 0 || i == sites - 1) {
      ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
          complex<double>(0.5, 0), MNud);
      ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
          complex<double>(0.5, 0), MNdu);
      ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
          complex<double>(Jz, 0), MNzz);
      ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
          complex<double>(-0.5 * Jz, 0), MNz);
      ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
          complex<double>(-0.5 * Jz, 0), MNzN);
      ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
          complex<double>(0.25 * Jz, 0), unitMn);
    }
    else {
    */
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(0.5, 0), MNud);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(-0.5, 0), MNdu);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(Jz, 0), MNzz);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(-0.5 * Jz, 0), MNz);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(-0.5 * Jz, 0), MNzN);
    ans += pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> >(
        complex<double>(0.25 * Jz, 0), unitMn);
  }
  return ans;
}

inline SpinHalfState1D Sud(size_t index, vector<bool> & Nums, double pref) {
  Nums[index] = true;
  Nums[index + 1] = false;
  return SpinHalfState1D(complex<double>(pref, 0), SpinHalfBaseState1D(Nums));
}

inline SpinHalfState1D Sdu(size_t index, vector<bool> & Nums, double pref) {
  Nums[index] = false;
  Nums[index + 1] = true;
  return SpinHalfState1D(complex<double>(pref, 0), SpinHalfBaseState1D(Nums));
}

inline SpinHalfState1D Szz(size_t index, vector<bool> & Nums, double pref) {
  if (Nums[index] == Nums[index + 1]) {
    return SpinHalfState1D(complex<double>(0.25 * pref, 0), SpinHalfBaseState1D(Nums));
  }
  else {
    return SpinHalfState1D(complex<double>(-0.25 * pref, 0), SpinHalfBaseState1D(Nums));
  }
}

SpinHalfState1D XXZ1D::makeMidState(size_t sites,
                                    double Jz,
                                    const SpinHalfBaseState1D & rhs) {
  SpinHalfState1D ans;
  vector<bool> Nums(rhs.getAllNums());
  for (size_t i = 0; i < sites - 1; i++) {
    if (i == 0 || i == sites - 2) {
      if (Nums[i] == false && Nums[i + 1] == true) {
        Nums[i] = true;
        Nums[i + 1] = false;
        ans += SpinHalfState1D(complex<double>(0.25, 0), SpinHalfBaseState1D(Nums));
        Nums[i] = false;
        Nums[i + 1] = true;
      }
      else if (Nums[i] == true && Nums[i + 1] == false) {
        Nums[i] = false;
        Nums[i + 1] = true;
        ans += SpinHalfState1D(complex<double>(0.25, 0), SpinHalfBaseState1D(Nums));
        Nums[i] = true;
        Nums[i + 1] = false;
      }
      if (Nums[i] == Nums[i + 1]) {
        ans += SpinHalfState1D(complex<double>(0.25 * 0.5 * Jz, 0),
                               SpinHalfBaseState1D(Nums));
      }
      else {
        ans += SpinHalfState1D(complex<double>(-0.25 * 0.5 * Jz, 0),
                               SpinHalfBaseState1D(Nums));
      }
      //ans += Szz(i, Nums, 0.5 * Jz);
    }
    else {
      if (Nums[i] == false && Nums[i + 1] == true) {
        Nums[i] = true;
        Nums[i + 1] = false;
        ans += SpinHalfState1D(complex<double>(0.5, 0), SpinHalfBaseState1D(Nums));
        Nums[i] = false;
        Nums[i + 1] = true;
      }
      else if (Nums[i] == true && Nums[i + 1] == false) {
        Nums[i] = false;
        Nums[i + 1] = true;
        ans += SpinHalfState1D(complex<double>(0.5, 0), SpinHalfBaseState1D(Nums));
        Nums[i] = true;
        Nums[i + 1] = false;
      }
      if (Nums[i] == Nums[i + 1]) {
        ans += SpinHalfState1D(complex<double>(0.25 * Jz, 0), SpinHalfBaseState1D(Nums));
      }
      else {
        ans += SpinHalfState1D(complex<double>(-0.25 * Jz, 0), SpinHalfBaseState1D(Nums));
      }
      //ans += Szz(i, Nums, Jz);
    }
  }
  return ans;
}

#endif  //ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
