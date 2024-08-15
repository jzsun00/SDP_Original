/*
  Jiazheng Sun
  Updated: Aug 8, 2024
  
  Class Implementations:
  Fermi1DOpSubBasis
  Fermi1DOpBasis
  
  Function Implementations:
  vector<pair<size_t, size_t> > findHermPairs(const Fermi1DOpBasis & basis)
  string printHermPairs(const vector<pair<size_t, size_t> > & pairs)
  void transVecToReIm(vector<complex<double> > & vec, vector<pair<size_t, size_t> > & pairs)
*/

#ifndef QM_FERMI_SUBSPACES_NONTEM_CPP
#define QM_FERMI_SUBSPACES_NONTEM_CPP

#include <exception>
#include <set>
#include <stdexcept>

#include "./fermiSubspaces.hpp"

using std::complex;
using std::pair;
using std::set;
using std::vector;

//----------------------------------------------------------Fermi1DOpSubBasis------------

Fermi1DOpSubBasis & Fermi1DOpSubBasis::operator=(const Fermi1DOpSubBasis & rhs) {
  if (this != &rhs) {
    this->start = rhs.start;
    this->end = rhs.end;
    this->order = rhs.order;
    this->Basis = rhs.Basis;
  }
  return *this;
}

void FermiMonomialsGenerator(vector<int> & current,
                             int start,
                             int n,
                             int maxValue,
                             vector<vector<int> > & result) {
  if (n == 0) {
    result.push_back(current);
    return;
  }
  for (int i = start; i <= maxValue; ++i) {
    current.push_back(i);
    FermiMonomialsGenerator(current, i, n - 1, maxValue, result);
    current.pop_back();
  }
}

FermiMonomial<Fermi1DLadderOp> FermiIntToMn(vector<int> input, bool creatorF) {
  FermiMonomial<Fermi1DLadderOp> ans;
  for (size_t index = 0; index < input.size(); index++) {
    Fermi1DLadderOp op(input[index], creatorF);
    ans *= op;
  }
  return ans;
}

vector<FermiMonomial<Fermi1DLadderOp> > FermiListMonomials(size_t length,
                                                           int start,
                                                           int end,
                                                           bool creatorF) {
  vector<int> current;
  vector<vector<int> > result;
  FermiMonomialsGenerator(current, start, length, end, result);
  vector<FermiMonomial<Fermi1DLadderOp> > ans;
  for (size_t i = 0; i < result.size(); ++i) {
    ans.push_back(FermiIntToMn(result[i], creatorF));
  }
  return ans;
}

void Fermi1DOpSubBasis::init(bool isInf) {
  for (size_t m = 0; m <= order; ++m) {
    if (m != order / 2) {
      continue;
    }
    if (m == 0) {
      vector<FermiMonomial<Fermi1DLadderOp> > creation =
          FermiListMonomials(order, start, end, true);
      for (size_t i = 0; i < creation.size(); i++) {
        creation[i].reverse();
      }
      Basis.insert(Basis.end(), creation.begin(), creation.end());
    }
    else if (m == order) {
      vector<FermiMonomial<Fermi1DLadderOp> > annihilation =
          FermiListMonomials(order, start, end, false);
      Basis.insert(Basis.end(), annihilation.begin(), annihilation.end());
    }
    else {
      vector<FermiMonomial<Fermi1DLadderOp> > creation =
          FermiListMonomials(order - m, start, end, true);
      for (size_t i = 0; i < creation.size(); i++) {
        creation[i].reverse();
      }
      vector<FermiMonomial<Fermi1DLadderOp> > annihilation =
          FermiListMonomials(m, start, end, false);
      for (size_t i = 0; i < annihilation.size(); ++i) {
        for (size_t j = 0; j < creation.size(); ++j) {
          FermiMonomial<Fermi1DLadderOp> copy(annihilation[i]);
          copy *= creation[j];
          if (isInf) {
            if (isNew(copy)) {
              Basis.push_back(copy);
            }
          }
          else {
            Basis.push_back(copy);
          }
        }
      }
    }
  }
}

std::string Fermi1DOpSubBasis::toString() const {
  std::string ans;
  ans += "Number of basis operators = ";
  ans += std::to_string(Basis.size());
  ans += "\nFull Basis:\n";
  size_t count = 1;
  for (auto it = Basis.begin(); it != Basis.end(); ++it) {
    ans += (std::to_string(count) + "\t");
    ans += it->toString();
    ans += "\n";
    count++;
  }
  return ans;
}

bool Fermi1DOpSubBasis::isNew(const FermiMonomial<Fermi1DLadderOp> & toAdd) const {
  const size_t len = Basis.size();
  for (size_t i = 0; i < len; i++) {
    if (toAdd.equiv(Basis[i])) {
      return false;
    }
  }
  return true;
}

//---------------------------------------------------------------Fermi1DOpBasis-------

Fermi1DOpBasis::Fermi1DOpBasis() : OpBasis<FermiMonomial<Fermi1DLadderOp>, int>() {
  Fermi1DLadderOp unit(true);
  Basis.push_back(unit);
}

Fermi1DOpBasis & Fermi1DOpBasis::operator=(const Fermi1DOpBasis & rhs) {
  if (this != &rhs) {
    this->Basis = rhs.Basis;
    this->lookupTable = rhs.lookupTable;
  }
  return *this;
}

void Fermi1DOpBasis::buildTable() {
  const size_t len = Basis.size();
  for (size_t i = 0; i < len; i++) {
    lookupTable[Basis[i]] = i;
  }
}

std::string Fermi1DOpBasis::toString() const {
  std::string ans;
  ans += "Number of basis operators = ";
  ans += std::to_string(Basis.size());
  ans += "\nFull Basis:\n";
  size_t count = 1;
  for (vector<FermiMonomial<Fermi1DLadderOp> >::const_iterator it = Basis.begin();
       it != Basis.end();
       ++it) {
    ans += (std::to_string(count) + "    ");
    ans += it->toString();
    ans += "\n";
    count++;
  }
  return ans;
}

void Fermi1DOpBasis::addSubspace(
    const OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int> & rhs) {
  vector<FermiMonomial<Fermi1DLadderOp> > sub = rhs.getFullBasis();
  Basis.insert(Basis.end(), sub.begin(), sub.end());
}

size_t Fermi1DOpBasis::findIndex(const FermiMonomial<Fermi1DLadderOp> & mn) const {
  auto it = lookupTable.find(mn);
  if (it == lookupTable.end()) {
    std::string error("ERROR: Input Fermi Monomial ");
    error += mn.toString();
    error += " does not exist in the Basis!\n";
    throw std::runtime_error(error);
  }
  return it->second;
}

vector<complex<double> > Fermi1DOpBasis::projPolyInf(
    const FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > & poly) const {
  vector<complex<double> > ans(Basis.size());
  for (size_t index = 0; index < poly.getSize(); index++) {
    ans[findIndex(poly[index].second)] = poly[index].first;
  }
  return ans;
}

std::vector<size_t> Fermi1DOpBasis::projPolyInf(
    std::vector<std::complex<double> > & vec,
    const FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > & poly) const {
  vector<size_t> validIdx;
  for (size_t polyIdx = 0; polyIdx < poly.getSize(); polyIdx++) {
    FermiMonomial<Fermi1DLadderOp> cp1(poly[polyIdx].second);
    FermiMonomial<Fermi1DLadderOp> cp2(poly[polyIdx].second);
    cp1.moveIndex(-cp1[0].getIndex());
    cp2.moveIndex(-cp2[cp2.getSize() - 1].getIndex());
    int basisIdx;
    bool validPoly = false;
    try {
      basisIdx = findIndex(cp1);
      validPoly = true;
    }
    catch (std::exception & e) {
    }
    try {
      basisIdx = findIndex(cp2);
      validPoly = true;
    }
    catch (std::exception & e) {
    }
    if (!validPoly) {
      std::string error("ERROR: Input Fermi Monomial");
      error += " does not exist in the Basis!\n";
      throw std::runtime_error(error);
    }
    vec[basisIdx] = poly[polyIdx].first;
    validIdx.push_back(basisIdx);
  }
  return validIdx;
}

std::vector<size_t> Fermi1DOpBasis::projPolyFinite(
    std::vector<std::complex<double> > & vec,
    const FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > & poly) const {
  vector<size_t> validIdx;
  for (size_t polyIdx = 0; polyIdx < poly.getSize(); polyIdx++) {
    int basisIdx = findIndex(poly[polyIdx].second);
    vec[basisIdx] = poly[polyIdx].first;
    validIdx.push_back(basisIdx);
  }
  return validIdx;
}

//--------------------------------------------------------------Other Functions----------

vector<pair<size_t, size_t> > FermiFindHermPairs(const Fermi1DOpBasis & basis) {
  set<size_t> addedIndex;  //Index already scanned or added
  vector<pair<size_t, size_t> > ans;
  const size_t len = basis.getLength();
  for (size_t index = 0; index < len; index++) {
    if (addedIndex.find(index) != addedIndex.end()) {  //Already scanned index
      continue;
    }
    addedIndex.insert(index);
    FermiMonomial<Fermi1DLadderOp> current = basis[index];
    FermiMonomial<Fermi1DLadderOp> currentCopy(current);
    currentCopy.herm();
    if (currentCopy == current) {  //current is Hermitian
      continue;
    }
    size_t hermIndex = basis.findIndex(currentCopy);  //Index of current.herm()
    addedIndex.insert(hermIndex);
    ans.push_back(pair<size_t, size_t>(index, hermIndex));
  }
  return ans;
}

std::string FermiPrintHermPairs(const vector<pair<size_t, size_t> > & pairs) {
  std::string ans = "";
  for (size_t i = 0; i < pairs.size(); i++) {
    ans += std::to_string(i + 1);
    ans += "    ";
    ans += std::to_string(pairs[i].first + 1);
    ans += ",  ";
    ans += std::to_string(pairs[i].second + 1);
    ans += "\n";
  }
  return ans;
}

void FermiTransVecToReIm(std::vector<std::complex<double> > & vec,
                         const std::vector<std::pair<size_t, size_t> > & pairs) {
  const size_t len = pairs.size();
  for (size_t i = 0; i < len; i++) {
    complex<double> ori1 = vec[pairs[i].first];
    complex<double> ori2 = vec[pairs[i].second];
    vec[pairs[i].first] = ori1 + ori2;
    vec[pairs[i].second] = complex<double>(0, 1.0) * (ori1 - ori2);
  }
}

#endif  //QM_FERMI_SUBSPACES_NONTEM_CPP
