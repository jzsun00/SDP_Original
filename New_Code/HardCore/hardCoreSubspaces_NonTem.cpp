/*
  Jiazheng Sun
  Updated: May 10, 2024

  Implementations of methods in class:
  Fermi1DLadderOp, FermiMonomial, FermiPolynomial.
 */

#ifndef ORI_SDP_GS_HARDCORESUBSPACES_NONTEM_CPP
#define ORI_SDP_GS_HARDCORESUBSPACES_NONTEM_CPP

#include <cstddef>

#include "./hardCoreSubspaces.hpp"

//-------------------------------------------------------------HardCore1DOpSubBasis------

void MonomialsGenerator(vector<int> & current,
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
    MonomialsGenerator(current, i, n - 1, maxValue, result);
    current.pop_back();
  }
}

HardCoreMonomial<HardCore1DLadderOp> intToMn(vector<int> input, bool creatorF) {
  HardCoreMonomial<HardCore1DLadderOp> ans;
  for (size_t index = 0; index < input.size(); index++) {
    HardCore1DLadderOp op(input[index], creatorF);
    ans *= op;
  }
  return ans;
}

vector<HardCoreMonomial<HardCore1DLadderOp> > listMonomials(size_t length,
                                                            int start,
                                                            int end,
                                                            bool creatorF) {
  vector<int> current;
  vector<vector<int> > result;
  MonomialsGenerator(current, start, length, end, result);
  vector<HardCoreMonomial<HardCore1DLadderOp> > ans;
  for (size_t i = 0; i < result.size(); ++i) {
    ans.push_back(intToMn(result[i], creatorF));
  }
  return ans;
}

void HardCore1DOpSubBasis::init() {
  for (size_t m = 0; m <= order; ++m) {
    if (m != order / 2) {
      continue;
    }
    if (m == 0) {
      vector<HardCoreMonomial<HardCore1DLadderOp> > creation =
          listMonomials(order, start, end, true);
      Basis.insert(Basis.end(), creation.begin(), creation.end());
    }
    else if (m == order) {
      vector<HardCoreMonomial<HardCore1DLadderOp> > annihilation =
          listMonomials(order, start, end, false);
      Basis.insert(Basis.end(), annihilation.begin(), annihilation.end());
    }
    else {
      vector<HardCoreMonomial<HardCore1DLadderOp> > creation =
          listMonomials(order - m, start, end, true);
      for (size_t i = 0; i < creation.size(); i++) {
        creation[i].reverse();
      }
      vector<HardCoreMonomial<HardCore1DLadderOp> > annihilation =
          listMonomials(m, start, end, false);
      for (size_t i = 0; i < annihilation.size(); ++i) {
        for (size_t j = 0; j < creation.size(); ++j) {
          HardCoreMonomial<HardCore1DLadderOp> copy(annihilation[i]);
          copy *= creation[j];
          //Basis.push_back(copy *= creation[j]);
          if (isNew(copy)) {
            Basis.push_back(copy);
          }
        }
      }
    }
  }
}

std::string HardCore1DOpSubBasis::toString() {
  std::string ans;
  ans += "Number of basis operators = ";
  ans += std::to_string(Basis.size());
  ans += "\nFull Basis:\n";
  size_t count = 1;
  for (vector<HardCoreMonomial<HardCore1DLadderOp> >::const_iterator it = Basis.begin();
       it != Basis.end();
       ++it) {
    ans += (std::to_string(count) + "    ");
    ans += it->toString();
    ans += "\n";
    count++;
  }
  return ans;
}

bool HardCore1DOpSubBasis::isNew(HardCoreMonomial<HardCore1DLadderOp> const & mn) {
  for (size_t i = 0; i < Basis.size(); i++) {
    if (mn.equiv(Basis[i])) {
      return false;
    }
  }
  return true;
}

//---------------------------------------------------------------HardCore1DOpBasis-------

std::string HardCore1DOpBasis::toString() {
  std::string ans;
  ans += "Number of basis operators = ";
  ans += std::to_string(Basis.size());
  ans += "\nFull Basis:\n";
  size_t count = 1;
  for (vector<HardCoreMonomial<HardCore1DLadderOp> >::const_iterator it = Basis.begin();
       it != Basis.end();
       ++it) {
    ans += (std::to_string(count) + "    ");
    ans += it->toString();
    ans += "\n";
    count++;
  }
  return ans;
}

void HardCore1DOpBasis::addSubspace(
    OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int> & rhs) {
  vector<HardCoreMonomial<HardCore1DLadderOp> > sub = rhs.getBasis();
  Basis.insert(Basis.end(), sub.begin(), sub.end());
}

vector<complex<double> > HardCore1DOpBasis::projPolyInf(
    HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly) {
  vector<complex<double> > ans(Basis.size());
  for (size_t index = 0; index < Basis.size(); index++) {
    HardCoreMonomial<HardCore1DLadderOp> basisMn = Basis[index];
    for (typename vector<pair<complex<double>, HardCoreMonomial<HardCore1DLadderOp> > >::
             const_iterator it = poly.getBegin();
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

vector<pair<size_t, size_t> > findHermPairs(HardCore1DOpBasis & basis) {
  set<size_t> addedIndex;
  vector<pair<size_t, size_t> > ans;
  for (size_t index = 1; index < basis.getLength(); index++) {
    if (addedIndex.find(index) != addedIndex.end()) {
      continue;
    }
    addedIndex.insert(index);
    HardCoreMonomial<HardCore1DLadderOp> current = basis[index];
    HardCoreMonomial<HardCore1DLadderOp> currentCopy(current);
    currentCopy.herm();
    if (currentCopy == current) {
      continue;
    }
    HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > PolyCurrent(current);
    HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > PolyCurrentCopy(
        currentCopy);
    PolyCurrentCopy.normalize();
    if (PolyCurrent == PolyCurrentCopy) {
      continue;
    }
    for (size_t j = index; j < basis.getLength(); j++) {
      HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > PolyBasis(basis[j]);
      if (PolyCurrentCopy == PolyBasis) {
        addedIndex.insert(j);
        ans.push_back(pair<size_t, size_t>(index, j));
      }
    }
  }
  return ans;
}

std::string printHermPairs(vector<pair<size_t, size_t> > & pairs) {
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

void transVecToReIm(vector<complex<double> > & vec,
                    vector<pair<size_t, size_t> > & pairs) {
  for (size_t i = 0; i < pairs.size(); i++) {
    complex<double> ori1 = vec[pairs[i].first];
    complex<double> ori2 = vec[pairs[i].second];
    vec[pairs[i].first] = ori1 + ori2;
    vec[pairs[i].second] = complex<double>(0, 1.0) * (ori1 - ori2);
  }
}

#endif  //ORI_SDP_GS_HARDCORESUBSPACES_NONTEM_CPP
