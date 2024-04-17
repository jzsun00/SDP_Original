/*
  Jiazheng Sun
  Updated: Apr 16, 2024

  Implementations of methods in class:
  Fermi1DLadderOp, FermiMonomial, FermiPolynomial.
 */

#ifndef ORI_SDP_GS_HARDCORECONSTRAINTS_NONTEM_CPP
#define ORI_SDP_GS_HARDCORECONSTRAINTS_NONTEM_CPP

#include <fstream>

#include "./hardCoreConstraints.hpp"

//-------------------------------------------------------------HardCore1DOpSubBasis------

void ConsMonomialsGenerator(vector<int> & current,
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
    ConsMonomialsGenerator(current, i, n - 1, maxValue, result);
    current.pop_back();
  }
}

HardCoreMonomial<HardCore1DLadderOp> ConsintToMn(vector<int> input, bool creatorF) {
  HardCoreMonomial<HardCore1DLadderOp> ans;
  for (size_t index = 0; index < input.size(); index++) {
    HardCore1DLadderOp op(input[index], creatorF);
    ans *= op;
  }
  return ans;
}

vector<HardCoreMonomial<HardCore1DLadderOp> > ConslistMonomials(size_t length,
                                                                int start,
                                                                int end,
                                                                bool creatorF) {
  vector<int> current;
  vector<vector<int> > result;
  ConsMonomialsGenerator(current, start, length, end, result);
  vector<HardCoreMonomial<HardCore1DLadderOp> > ans;
  for (size_t i = 0; i < result.size(); ++i) {
    ans.push_back(ConsintToMn(result[i], creatorF));
  }
  return ans;
}

void HardCore1DConsBaseSet::init() {
  for (size_t m = 0; m <= order; ++m) {
    if (m == 0) {
      continue;
    }
    if (m == 0) {
      vector<HardCoreMonomial<HardCore1DLadderOp> > creation =
          ConslistMonomials(order, start, end, true);
      BaseOpSet.insert(BaseOpSet.end(), creation.begin(), creation.end());
    }
    else if (m == order) {
      vector<HardCoreMonomial<HardCore1DLadderOp> > annihilation =
          ConslistMonomials(order, start, end, false);
      BaseOpSet.insert(BaseOpSet.end(), annihilation.begin(), annihilation.end());
    }
    else {
      vector<HardCoreMonomial<HardCore1DLadderOp> > creation =
          ConslistMonomials(order - m, start, end, true);
      vector<HardCoreMonomial<HardCore1DLadderOp> > annihilation =
          ConslistMonomials(m, start, end, false);
      for (size_t i = 0; i < annihilation.size(); ++i) {
        for (size_t j = 0; j < creation.size(); ++j) {
          HardCoreMonomial<HardCore1DLadderOp> copy(annihilation[i]);
          BaseOpSet.push_back(copy *= creation[j]);
        }
      }
    }
  }
}

std::string HardCore1DConsBaseSet::toString() {
  std::string ans;
  ans += "Number of basis operators = ";
  ans += std::to_string(BaseOpSet.size());
  ans += "\nFull Basis:\n";
  size_t count = 1;
  for (vector<HardCoreMonomial<HardCore1DLadderOp> >::const_iterator it =
           BaseOpSet.begin();
       it != BaseOpSet.end();
       ++it) {
    ans += (std::to_string(count) + "    ");
    ans += it->toString();
    ans += "\n";
    count++;
  }
  return ans;
}

//------------------------------------------------------------HardCore1DConsSet----------

std::string HardCore1DConsSet::toString() {
  std::string ans;
  ans += "Number of basis operators = ";
  ans += std::to_string(OpSet.size());
  ans += "\nFull Basis:\n";
  size_t count = 1;
  for (vector<HardCoreMonomial<HardCore1DLadderOp> >::const_iterator it = OpSet.begin();
       it != OpSet.end();
       ++it) {
    ans += (std::to_string(count) + "    ");
    ans += it->toString();
    ans += "\n";
    count++;
  }
  return ans;
}

void HardCore1DConsSet::addBaseSet(
    ConsBaseSet<HardCoreMonomial<HardCore1DLadderOp>, int> & rhs) {
  vector<HardCoreMonomial<HardCore1DLadderOp> > sub = rhs.getBaseOpSet();
  OpSet.insert(OpSet.end(), sub.begin(), sub.end());
}

HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > HardCore1DConsSet::getIJPoly(
    size_t i,
    size_t j) {
  HardCoreMonomial<HardCore1DLadderOp> mnI = OpSet[i];
  mnI.herm();
  HardCoreMonomial<HardCore1DLadderOp> mnJ = OpSet[j];
  mnI *= mnJ;
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > ans(mnI);
  ans.normalize();
  return ans;
}

//-------------------------------------------------------------Other Functions-----------

void printMatrixHardCore1D(HardCore1DConsSet & constraints, HardCore1DOpBasis & basis, std::string fileName, vector<complex<double> > ham) {
  size_t matrixNum = basis.getLength();
  size_t matrixSize = constraints.getLength();
  vector<vector<vector<complex<double> > > > matrices(
      matrixNum,
      vector<vector<complex<double> > >(matrixSize,
                                        vector<complex<double> >(matrixSize)));
  for (size_t i = 0; i < matrixSize; i++) {
    for (size_t j = 0; j < matrixSize; j++) {
      HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > polyIJ =
          constraints.getIJPoly(i, j);
      vector<complex<double> > entryIJ = basis.projPoly(polyIJ);
      //std::cout << "i = " << i << ",  j = " << j << std::endl;
      //std::cout << "polyIJ = " << polyIJ.toString() << std::endl;
      //std::cout << "entryIJ = " << complexVector_toString(entryIJ) << std::endl;
      for (size_t k = 0; k < matrixNum; k++) {
        matrices[k][i][j] = entryIJ[k];
      }
    }
  }
  for (size_t num = 0; num < matrixNum; num++) {
    //std::cout << "num = " << num << std::endl;
    //std::cout << basis[num].toString() << std::endl;
    //std::cout << complexMatrix_toString(matrices[num]) << std::endl;
  }
  std::ofstream inputFile(fileName);
  if (!inputFile.is_open()) {
    std::cerr << "Failed to open file for writing." << std::endl;
  }
  inputFile << "\"XXZ Test 1: mDim = 16, nBLOCK = 1, {6}\"" << std::endl;
  inputFile << matrixNum - 1 <<"  =  mDIM" << std::endl;
  inputFile << "1  =  nBLOCK" << std::endl;
  inputFile << matrixSize - 1 << "  = bLOCKsTRUCT" << std::endl;
  inputFile << "{";
  for (size_t i = 1; i < ham.size(); i++) {
    inputFile << ham[i].real();
    if (i < ham.size() - 1) {
      inputFile << ", ";
    }
  }
  inputFile << " }" << std::endl;
  for (size_t num = 0; num < matrixNum; num++) {
    inputFile << "{";
    for (size_t i = 1; i < matrixSize; i++) {
      inputFile << "{";
      for (size_t j = 1; j < matrixSize; j++) {
        if (num == 0) {
          inputFile << -1 * matrices[num][i][j].real();
        }
        else {
          inputFile << matrices[num][i][j].real();
        }
        if (j < matrixSize - 1) {
          inputFile << ", ";
        }
      }
      inputFile << "} ";
    }
    inputFile << "}" << std::endl;
  }
  inputFile.close();
  std::cout << "File has been written successfully." << std::endl;
}

#endif  //ORI_SDP_GS_HARDCORECONSTRAINTS_NONTEM_CPP
