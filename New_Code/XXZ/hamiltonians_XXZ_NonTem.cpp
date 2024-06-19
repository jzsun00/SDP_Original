/*
  Jiazheng Sun
  Updated: Jun 18, 2024

  Define dense and sparse Hamiltonian matrices for 1D XXZ model.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP
#define ORI_SDP_GS_HAMILTONIANS_XXZ_NONTEM_CPP

#include <omp.h>

#include <chrono>
#include <cstddef>
#include <map>
#include <vector>

#include "./hamiltonians_XXZ.hpp"

//------------------------------------------------------------XXZSparseHamiltonian-------

void XXZSparseHamiltonian::createMatrix(SpinHalfBasis1D & basis) {
  dim = basis.getSize();
  pcol.push_back(0);
  vector<SpinHalfState1D> midStates(dim);

  auto start_calc_all_mid = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
  for (long unsigned j = 0; j < dim; j++) {
    midStates[j] = makeMidState(sites, Jz, basis[j]);
  }
  auto end_calc_all_mid = std::chrono::high_resolution_clock::now();
  auto duration_calc_all_mid = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_calc_all_mid - start_calc_all_mid);
  std::cout << "\nCalculating all mid Running Time: " << duration_calc_all_mid.count()
            << " ms" << std::endl;

  //Use j = 0 to test performance
  for (long unsigned j = 0; j < 2; j++) {
    auto start_calc_mid = std::chrono::high_resolution_clock::now();
    //SpinHalfState1D mid = poly * basis[j];
    //SpinHalfState1D mid = makeMidState(sites, Jz, basis[j]);
    //SpinHalfState1D mid = midStates[j];
    auto end_calc_mid = std::chrono::high_resolution_clock::now();
    //mid.eraseZeros();
    std::vector<int> indices;
    std::map<int, complex<double> > elements;
    auto start_fill_map = std::chrono::high_resolution_clock::now();
    for (size_t k = 0; k < midStates[j].getSize(); k++) {
      //int index = basis.findBaseState(mid[k].second);
      size_t index = basis.lookUpBaseState(midStates[j][k].second);
      complex<double> value = midStates[j][k].first;
      indices.push_back(index);
      elements[index] = value;
    }
    auto end_fill_map = std::chrono::high_resolution_clock::now();
    std::sort(indices.begin(), indices.end());
    auto start_fill_matrix = std::chrono::high_resolution_clock::now();
    for (size_t k = 0; k < midStates[j].getSize(); k++) {
      nnz++;
      nzVal.push_back(elements[indices[k]]);
      irow.push_back(indices[k]);
    }
    auto end_fill_matrix = std::chrono::high_resolution_clock::now();
    pcol.push_back(nnz);
    // Record time.
    auto duration_calc_mid = std::chrono::duration_cast<std::chrono::nanoseconds>(
        end_calc_mid - start_calc_mid);
    auto duration_fill_map = std::chrono::duration_cast<std::chrono::nanoseconds>(
        end_fill_map - start_fill_map);
    auto duration_fill_matrix = std::chrono::duration_cast<std::chrono::nanoseconds>(
        end_fill_matrix - start_fill_matrix);
    auto duration_total = std::chrono::duration_cast<std::chrono::nanoseconds>(
        end_fill_matrix - start_calc_mid);
    std::cout << "\nCalculating mid Running Time: " << duration_calc_mid.count() << " ns"
              << std::endl;
    std::cout << "\nFilling Map Running Time: " << duration_fill_map.count() << " ns"
              << std::endl;
    std::cout << "\nFilling Matrix Running Time: " << duration_fill_matrix.count()
              << " ns" << std::endl;
    std::cout << "\nTotal Running Time: " << duration_total.count() << " ns" << std::endl;
  }

  //Continue from j = 1
  for (long unsigned j = 2; j < dim; j++) {
    //SpinHalfState1D mid = poly * basis[j];
    //SpinHalfState1D mid = makeMidState(sites, Jz, basis[j]);
    //SpinHalfState1D mid = midStates[j];
    //mid.eraseZeros();
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

SpinHalfState1D makeMidState(size_t sites, double Jz, const SpinHalfBaseState1D & rhs) {
  SpinHalfState1D ans;
  vector<bool> Nums(rhs.getNums());
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
