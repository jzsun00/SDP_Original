/*
  Jiazheng Sun
  Updated: Aug 6, 2024
  
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

#include <set>

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

void Fermi1DOpSubBasis::init() {
  for (size_t m = 0; m <= order; ++m) {
    if (m != order / 2) {
      continue;
    }
    if (m == 0) {
      vector<FermiMonomial<Fermi1DLadderOp> > creation =
          FermiListMonomials(order, start, end, true);
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
          //if (isNew(copy)) {
          Basis.push_back(copy);
          //}
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
  size_t len = Basis.size();
  for (size_t i = 0; i < len; i++) {
    if (toAdd.equiv(Basis[i])) {
      return false;
    }
  }
  return true;
}

//---------------------------------------------------------------Fermi1DOpBasis-------

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

vector<complex<double> > Fermi1DOpBasis::projPolyInf(
    FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly) {
  vector<complex<double> > ans(Basis.size());
  for (size_t index = 0; index < Basis.size(); index++) {
    FermiMonomial<Fermi1DLadderOp> basisMn = Basis[index];
    for (typename vector<
             pair<complex<double>, FermiMonomial<Fermi1DLadderOp> > >::const_iterator it =
             poly.getBegin();
         it != poly.getEnd();
         ++it) {
      if (it->second.equiv(basisMn)) {
        ans[index] = it->first;
      }
    }
  }
  return ans;
}

//--------------------------------------------------------------Other Functions----------

vector<pair<size_t, size_t> > FermiFindHermPairs(Fermi1DOpBasis & basis) {
  set<size_t> addedIndex;
  vector<pair<size_t, size_t> > ans;
  for (size_t index = 1; index < basis.getLength(); index++) {
    if (addedIndex.find(index) != addedIndex.end()) {
      continue;
    }
    addedIndex.insert(index);
    FermiMonomial<Fermi1DLadderOp> current = basis[index];
    FermiMonomial<Fermi1DLadderOp> currentCopy(current);
    currentCopy.herm();
    if (currentCopy == current) {
      continue;
    }
    FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > PolyCurrent(current);
    FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > PolyCurrentCopy(currentCopy);
    PolyCurrentCopy.normalOrder();
    if (PolyCurrent == PolyCurrentCopy) {
      continue;
    }
    for (size_t j = index; j < basis.getLength(); j++) {
      FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > PolyBasis(basis[j]);
      if (PolyCurrentCopy == PolyBasis) {
        addedIndex.insert(j);
        ans.push_back(pair<size_t, size_t>(index, j));
      }
    }
  }
  return ans;
}

std::string FermiPrintHermPairs(vector<pair<size_t, size_t> > & pairs) {
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
                         std::vector<std::pair<size_t, size_t> > & pairs) {
  for (size_t i = 0; i < pairs.size(); i++) {
    complex<double> ori1 = vec[pairs[i].first];
    complex<double> ori2 = vec[pairs[i].second];
    vec[pairs[i].first] = ori1 + ori2;
    vec[pairs[i].second] = complex<double>(0, 1.0) * (ori1 - ori2);
  }
}

#endif  //QM_FERMI_SUBSPACES_NONTEM_CPP
