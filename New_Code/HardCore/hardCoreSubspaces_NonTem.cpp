/*
  Jiazheng Sun
  Updated: Apr 6, 2024

  Implementations of methods in class:
  Fermi1DLadderOp, FermiMonomial, FermiPolynomial.
 */

#ifndef ORI_SDP_GS_HARDCORESUBSPACES_NONTEM_CPP
#define ORI_SDP_GS_HARDCORESUBSPACES_NONTEM_CPP

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
      vector<HardCoreMonomial<HardCore1DLadderOp> > annihilation =
          listMonomials(m, start, end, false);
      for (size_t i = 0; i < annihilation.size(); ++i) {
        for (size_t j = 0; j < creation.size(); ++j) {
          HardCoreMonomial<HardCore1DLadderOp> copy(annihilation[i]);
          Basis.push_back(copy *= creation[j]);
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

#endif  //ORI_SDP_GS_HARDCORESUBSPACES_NONTEM_CPP
