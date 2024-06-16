/*
  Jiazheng Sun
  Updated: Apr 16, 2024

  Implementations of methods in class:
  Fermi1DLadderOp, FermiMonomial, FermiPolynomial.
 */

#ifndef ORI_SDP_GS_HARDCORECONSTRAINTS_NONTEM_CPP
#define ORI_SDP_GS_HARDCORECONSTRAINTS_NONTEM_CPP

#include <fstream>
#include <iostream>

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
    //if (m == 0) {
    //  continue;
    //}
    if (m == 0) {
      vector<HardCoreMonomial<HardCore1DLadderOp> > creation =
          ConslistMonomials(order, start, end, true);
      for (size_t i = 0; i < creation.size(); i++) {
        creation[i].reverse();
      }
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
  ans += "\nFull Constraint Set:\n";
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

void printMatrixHardCore1D(HardCore1DConsSet & constraints,
                           HardCore1DOpBasis & basis,
                           std::string fileName,
                           vector<complex<double> > ham,
                           vector<pair<size_t, size_t> > & pairs) {
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
  transMatToReIm(matrices, pairs);
  std::ofstream inputFile(fileName);
  if (!inputFile.is_open()) {
    std::cerr << "Failed to open file for writing." << std::endl;
  }
  inputFile << "\"XXZ Test: mDim = " << matrixNum - 1 << ", nBLOCK = 1, {"
            << matrixSize * 2 << "}\"" << std::endl;
  inputFile << matrixNum - 1 << "  =  mDIM" << std::endl;
  inputFile << "1  =  nBLOCK" << std::endl;
  inputFile << matrixSize * 2 << "  = bLOCKsTRUCT" << std::endl;
  inputFile << "{";
  for (size_t i = 1; i < ham.size(); i++) {
    inputFile << ham[i].real();
    if (i < ham.size() - 1) {
      inputFile << ", ";
    }
  }
  inputFile << " }" << std::endl;
  for (size_t num = 0; num < matrixNum; num++) {
    inputFile << "{ ";
    for (size_t i = 0; i < matrixSize; i++) {
      inputFile << "{";
      // First Block
      for (size_t j = 0; j < matrixSize; j++) {
        if (num == 0) {
          inputFile << -matrices[num][i][j].real();
        }
        else {
          inputFile << matrices[num][i][j].real();
        }
        inputFile << ", ";
      }
      // Second Block
      for (size_t j = 0; j < matrixSize; j++) {
        if (num == 0) {
          inputFile << matrices[num][i][j].imag();
        }
        else {
          inputFile << -matrices[num][i][j].imag();
        }
        if (j < matrixSize - 1) {
          inputFile << ", ";
        }
      }
      inputFile << "},\n";
    }
    for (size_t i = 0; i < matrixSize; i++) {
      inputFile << "{";
      // Third Block
      for (size_t j = 0; j < matrixSize; j++) {
        if (num == 0) {
          inputFile << -matrices[num][i][j].imag();
        }
        else {
          inputFile << matrices[num][i][j].imag();
        }
        inputFile << ", ";
      }
      // Fourth Block
      for (size_t j = 0; j < matrixSize; j++) {
        if (num == 0) {
          inputFile << -matrices[num][i][j].real();
        }
        else {
          inputFile << matrices[num][i][j].real();
        }
        if (j < matrixSize - 1) {
          inputFile << ", ";
        }
      }
      if (i < matrixSize - 1) {
        inputFile << "},\n";
      }
      else {
        inputFile << "}";
      }
    }
    inputFile << " }" << std::endl;
  }
  inputFile.close();
  std::cout << "File has been written successfully." << std::endl;
}
////////////////////////////////////////////////////////////////////////////////////
void printMatrixXX1D(size_t max, std::string fileName) {
  size_t matrixNum = max * 2;
  size_t matrixSize = max * 2;
  vector<vector<vector<complex<double> > > > matrices(
      matrixNum,
      vector<vector<complex<double> > >(matrixSize,
                                        vector<complex<double> >(matrixSize)));
  for (size_t i = max; i < matrixSize; i++) {
    matrices[0][i][i] = -1;
  }
  for (size_t i = 0; i < max; i++) {
    matrices[1][i][i] = 1;
    matrices[1][i + max][i + max] = -1;
  }
  for (size_t dist = 1; dist < max; dist++) {
    size_t Num = max - dist;
    //std::cout << "dist = " << dist << ",  Num = " << Num << std::endl;
    for (size_t i = 0; i < Num; i++) {
      //std::cout << "i = " << i << std::endl;
      matrices[dist + 1][i][dist + i] = 1;
      matrices[dist + 1][dist + i][i] = 1;
      matrices[dist + 1][i + max][dist + i + max] = -1;
      matrices[dist + 1][dist + i + max][i + max] = -1;
      matrices[dist + max][i][dist + i] = complex<double>(0, -1);
      matrices[dist + max][dist + i][i] = complex<double>(0, 1);
      matrices[dist + max][i + max][dist + i + max] = complex<double>(0, -1);
      matrices[dist + max][dist + i + max][i + max] = complex<double>(0, 1);
    }
  }
  std::ofstream inputFile(fileName);
  if (!inputFile.is_open()) {
    std::cerr << "Failed to open file for writing." << std::endl;
  }
  inputFile << "\"XXZ Test: mDim = " << matrixNum - 1 << ", nBLOCK = 1, {"
            << matrixSize * 2 << "}\"" << std::endl;
  inputFile << matrixNum - 1 << "  =  mDIM" << std::endl;
  inputFile << "1  =  nBLOCK" << std::endl;
  inputFile << matrixSize * 2 << "  = bLOCKsTRUCT" << std::endl;
  inputFile << "{";
  for (size_t i = 1; i < 2 * max; i++) {
    if (i == 2) {
      inputFile << max - 1;
    }
    else {
      inputFile << 0;
    }
    if (i < 2 * max - 1) {
      inputFile << ", ";
    }
  }
  inputFile << " }" << std::endl;
  for (size_t num = 0; num < matrixNum; num++) {
    inputFile << "{ ";
    for (size_t i = 0; i < matrixSize; i++) {
      inputFile << "{";
      // First Block
      for (size_t j = 0; j < matrixSize; j++) {
        inputFile << matrices[num][i][j].real();
        inputFile << ", ";
      }
      // Second Block
      for (size_t j = 0; j < matrixSize; j++) {
        inputFile << -matrices[num][i][j].imag();
        if (j < matrixSize - 1) {
          inputFile << ", ";
        }
      }
      inputFile << "},\n";
    }
    for (size_t i = 0; i < matrixSize; i++) {
      inputFile << "{";
      // Third Block
      for (size_t j = 0; j < matrixSize; j++) {
        inputFile << matrices[num][i][j].imag();
        inputFile << ", ";
      }
      // Fourth Block
      for (size_t j = 0; j < matrixSize; j++) {
        inputFile << matrices[num][i][j].real();
        if (j < matrixSize - 1) {
          inputFile << ", ";
        }
      }
      if (i < matrixSize - 1) {
        inputFile << "},\n";
      }
      else {
        inputFile << "}";
      }
    }
    inputFile << " }" << std::endl;
  }
  inputFile.close();
  std::cout << "File has been written successfully." << std::endl;
}
void printSparseMatrixXX1D(size_t max, std::string fileName) {
  size_t matrixNum = max * 2;
  size_t matrixSize = max * 2;
  vector<vector<vector<complex<double> > > > matrices(
      matrixNum,
      vector<vector<complex<double> > >(matrixSize,
                                        vector<complex<double> >(matrixSize)));
  for (size_t i = max; i < matrixSize; i++) {
    matrices[0][i][i] = -1;
  }
  for (size_t i = 0; i < max; i++) {
    matrices[1][i][i] = 1;
    matrices[1][i + max][i + max] = -1;
  }
  for (size_t dist = 1; dist < max; dist++) {
    size_t Num = max - dist;
    //std::cout << "dist = " << dist << ",  Num = " << Num << std::endl;
    for (size_t i = 0; i < Num; i++) {
      //std::cout << "i = " << i << std::endl;
      matrices[dist + 1][i][dist + i] = 1;
      matrices[dist + 1][dist + i][i] = 1;
      matrices[dist + 1][i + max][dist + i + max] = -1;
      matrices[dist + 1][dist + i + max][i + max] = -1;
      matrices[dist + max][i][dist + i] = complex<double>(0, -1);
      matrices[dist + max][dist + i][i] = complex<double>(0, 1);
      matrices[dist + max][i + max][dist + i + max] = complex<double>(0, -1);
      matrices[dist + max][dist + i + max][i + max] = complex<double>(0, 1);
    }
  }
  std::ofstream inputFile(fileName);
  if (!inputFile.is_open()) {
    std::cerr << "Failed to open file for writing." << std::endl;
  }
  inputFile << "\"XXZ Test: mDim = " << matrixNum - 1 << ", nBLOCK = 1, {"
            << matrixSize * 2 << "}\"" << std::endl;
  inputFile << matrixNum - 1 << "  =  mDIM" << std::endl;
  inputFile << "1  =  nBLOCK" << std::endl;
  inputFile << matrixSize * 2 << "  = bLOCKsTRUCT" << std::endl;
  inputFile << "{";
  for (size_t i = 1; i < 2 * max; i++) {
    if (i == 2) {
      //inputFile << max - 1;
      inputFile << 1;
    }
    else {
      inputFile << 0;
    }
    if (i < 2 * max - 1) {
      inputFile << ", ";
    }
  }
  inputFile << " }" << std::endl;
  for (size_t num = 0; num < matrixNum; num++) {
    for (size_t i = 0; i < matrixSize; i++) {
      // First Block
      for (size_t j = i; j < matrixSize; j++) {
        if (std::abs(matrices[num][i][j].real()) < ERROR) {
          continue;
        }
        inputFile << num << " "
                  << "1 " << i + 1 << " " << j + 1 << " " << matrices[num][i][j].real()
                  << std::endl;
      }
      // Second Block
      for (size_t j = i; j < matrixSize; j++) {
        if (std::abs(matrices[num][i][j].imag()) < ERROR) {
          continue;
        }

        inputFile << num << " "
                  << "1 " << i + 1 << " " << j + 1 + matrixSize << " "
                  << -matrices[num][i][j].imag() << std::endl;
      }
    }
    for (size_t i = 0; i < matrixSize; i++) {
      // Fourth Block
      for (size_t j = i; j < matrixSize; j++) {
        if (std::abs(matrices[num][i][j].real()) < ERROR) {
          continue;
        }

        inputFile << num << " "
                  << "1 " << i + 1 + matrixSize << " " << j + 1 + matrixSize << " "
                  << matrices[num][i][j].real() << std::endl;
      }
    }
  }
  inputFile.close();
  std::cout << "File has been written successfully." << std::endl;
}
//////////////////////////////////////////////////////////////////////////////////
void printSparseMatrixHardCore1D(HardCore1DConsSet & constraints,
                                 HardCore1DOpBasis & basis,
                                 std::string fileName,
                                 vector<complex<double> > ham,
                                 vector<pair<size_t, size_t> > & pairs) {
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
      //if ((polyIJ.getSize() % 2) != 0) {
      //  continue;
      //}
      //vector<complex<double> > entryIJ = basis.projPoly(polyIJ);
      vector<complex<double> > entryIJ = basis.projPolyInf(polyIJ);
      for (size_t k = 0; k < matrixNum; k++) {
        matrices[k][i][j] = entryIJ[k];
      }
    }
  }
  transMatToReIm(matrices, pairs);
  std::cout << "\nMatrix construction completed" << std::endl
            << "Start writing file" << std::endl;
  std::ofstream inputFile(fileName);
  if (!inputFile.is_open()) {
    std::cerr << "Failed to open file for writing." << std::endl;
  }
  inputFile << "\"XXZ Test: mDim = " << matrixNum - 1 << ", nBLOCK = 1, {"
            << matrixSize * 2 << "}\"" << std::endl;
  inputFile << matrixNum - 1 << "  =  mDIM" << std::endl;
  inputFile << "1  =  nBLOCK" << std::endl;
  inputFile << matrixSize * 2 << "  = bLOCKsTRUCT" << std::endl;
  inputFile << "{";
  for (size_t i = 1; i < ham.size(); i++) {
    if (i == 1) {
      inputFile << 2 * ham[i].real();
    }
    else {
      inputFile << ham[i].real();
    }
    if (i < ham.size() - 1) {
      inputFile << ", ";
    }
  }
  inputFile << " }" << std::endl;
  for (size_t num = 0; num < matrixNum; num++) {
    for (size_t i = 0; i < matrixSize; i++) {
      // First Block
      for (size_t j = i; j < matrixSize; j++) {
        if (std::abs(matrices[num][i][j].real()) < ERROR) {
          continue;
        }
        if (num == 0) {
          inputFile << num << " "
                    << "1 " << i + 1 << " " << j + 1 << " " << -matrices[num][i][j].real()
                    << std::endl;
        }
        else {
          inputFile << num << " "
                    << "1 " << i + 1 << " " << j + 1 << " " << matrices[num][i][j].real()
                    << std::endl;
        }
      }
      // Second Block
      for (size_t j = i; j < matrixSize; j++) {
        if (std::abs(matrices[num][i][j].imag()) < ERROR) {
          continue;
        }
        if (num == 0) {
          inputFile << num << " "
                    << "1 " << i + 1 << " " << j + 1 + matrixSize << " "
                    << matrices[num][i][j].imag() << std::endl;
        }
        else {
          inputFile << num << " "
                    << "1 " << i + 1 << " " << j + 1 + matrixSize << " "
                    << -matrices[num][i][j].imag() << std::endl;
        }
      }
    }
    for (size_t i = 0; i < matrixSize; i++) {
      // Fourth Block
      for (size_t j = i; j < matrixSize; j++) {
        if (std::abs(matrices[num][i][j].real()) < ERROR) {
          continue;
        }
        if (num == 0) {
          inputFile << num << " "
                    << "1 " << i + 1 + matrixSize << " " << j + 1 + matrixSize << " "
                    << -matrices[num][i][j].real() << std::endl;
        }
        else {
          inputFile << num << " "
                    << "1 " << i + 1 + matrixSize << " " << j + 1 + matrixSize << " "
                    << matrices[num][i][j].real() << std::endl;
        }
      }
    }
  }
  inputFile.close();
  std::cout << "File has been written successfully." << std::endl;
}

void transMatToReIm(vector<vector<vector<complex<double> > > > & matrices,
                    vector<pair<size_t, size_t> > & pairs) {
  for (size_t n = 0; n < pairs.size(); n++) {
    for (size_t i = 0; i < matrices[0].size(); i++) {
      for (size_t j = 0; j < matrices[0].size(); j++) {
        complex<double> ori1 = matrices[pairs[n].first][i][j];
        complex<double> ori2 = matrices[pairs[n].second][i][j];
        matrices[pairs[n].first][i][j] = ori1 + ori2;
        matrices[pairs[n].second][i][j] = complex<double>(0, 1.0) * (ori1 - ori2);
      }
    }
  }
}

#endif  //ORI_SDP_GS_HARDCORECONSTRAINTS_NONTEM_CPP
