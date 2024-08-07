/*
  Jiazheng Sun
  Updated: Aug 7, 2024
  
  Class Implementations:
  Fermi1DConsBaseSet
  Fermi1DConsSet
  
  Function Implementations:
*/

#ifndef QM_FERMI_CONSTRAINTS_NONTEM_CPP
#define QM_FERMI_CONSTRAINTS_NONTEM_CPP

#include <exception>
#include <fstream>
#include <iostream>

#include "./fermiConstraints.hpp"

using std::complex;
using std::pair;
using std::vector;

//-------------------------------------------------------------Fermi1DOpSubBasis------

Fermi1DConsBaseSet & Fermi1DConsBaseSet::operator=(const Fermi1DConsBaseSet & rhs) {
  if (this != &rhs) {
    this->start = rhs.start;
    this->end = rhs.end;
    this->order = rhs.order;
    this->BaseOpSet = rhs.BaseOpSet;
  }
  return *this;
}

void FermiConsMonomialsGenerator(vector<int> & current,
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
    FermiConsMonomialsGenerator(current, i, n - 1, maxValue, result);
    current.pop_back();
  }
}

FermiMonomial<Fermi1DLadderOp> FermiConsintToMn(vector<int> input, bool creatorF) {
  FermiMonomial<Fermi1DLadderOp> ans;
  for (size_t index = 0; index < input.size(); index++) {
    Fermi1DLadderOp op(input[index], creatorF);
    ans *= op;
  }
  return ans;
}

vector<FermiMonomial<Fermi1DLadderOp> > FermiConslistMonomials(size_t length,
                                                               int start,
                                                               int end,
                                                               bool creatorF) {
  vector<int> current;
  vector<vector<int> > result;
  FermiConsMonomialsGenerator(current, start, length, end, result);
  vector<FermiMonomial<Fermi1DLadderOp> > ans;
  for (size_t i = 0; i < result.size(); ++i) {
    ans.push_back(FermiConsintToMn(result[i], creatorF));
  }
  return ans;
}

void Fermi1DConsBaseSet::init() {
  for (size_t m = 0; m <= order; ++m) {
    //if (m == 0) {
    //  continue;
    //}
    if (m == 0) {
      vector<FermiMonomial<Fermi1DLadderOp> > creation =
          FermiConslistMonomials(order, start, end, true);
      for (size_t i = 0; i < creation.size(); i++) {
        creation[i].reverse();
      }
      BaseOpSet.insert(BaseOpSet.end(), creation.begin(), creation.end());
    }
    else if (m == order) {
      vector<FermiMonomial<Fermi1DLadderOp> > annihilation =
          FermiConslistMonomials(order, start, end, false);
      BaseOpSet.insert(BaseOpSet.end(), annihilation.begin(), annihilation.end());
    }
    else {
      vector<FermiMonomial<Fermi1DLadderOp> > creation =
          FermiConslistMonomials(order - m, start, end, true);
      vector<FermiMonomial<Fermi1DLadderOp> > annihilation =
          FermiConslistMonomials(m, start, end, false);
      for (size_t i = 0; i < annihilation.size(); ++i) {
        for (size_t j = 0; j < creation.size(); ++j) {
          FermiMonomial<Fermi1DLadderOp> copy(annihilation[i]);
          copy *= creation[j];
          BaseOpSet.push_back(copy);
        }
      }
    }
  }
}

std::string Fermi1DConsBaseSet::toString() const {
  std::string ans;
  ans += "Number of basis operators = ";
  ans += std::to_string(BaseOpSet.size());
  ans += "\nFull Basis:\n";
  size_t count = 1;
  for (vector<FermiMonomial<Fermi1DLadderOp> >::const_iterator it = BaseOpSet.begin();
       it != BaseOpSet.end();
       ++it) {
    ans += (std::to_string(count) + "    ");
    ans += it->toString();
    ans += "\n";
    count++;
  }
  return ans;
}

//------------------------------------------------------------Fermi1DConsSet----------

Fermi1DConsSet::Fermi1DConsSet() : ConsSet<FermiMonomial<Fermi1DLadderOp>, int>() {
  Fermi1DLadderOp unit(true);
  OpSet.push_back(unit);
}

std::string Fermi1DConsSet::toString() const {
  std::string ans;
  ans += "Number of basis operators = ";
  ans += std::to_string(OpSet.size());
  ans += "\nFull Constraint Set:\n";
  size_t count = 1;
  for (vector<FermiMonomial<Fermi1DLadderOp> >::const_iterator it = OpSet.begin();
       it != OpSet.end();
       ++it) {
    ans += (std::to_string(count) + "    ");
    ans += it->toString();
    ans += "\n";
    count++;
  }
  return ans;
}

Fermi1DConsSet & Fermi1DConsSet::operator=(const Fermi1DConsSet & rhs) {
  if (this != &rhs) {
    this->OpSet = rhs.OpSet;
  }
  return *this;
}

void Fermi1DConsSet::addBaseSet(
    const ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int> & rhs) {
  vector<FermiMonomial<Fermi1DLadderOp> > sub = rhs.getFullBaseOpSet();
  OpSet.insert(OpSet.end(), sub.begin(), sub.end());
}

FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > Fermi1DConsSet::getIJPoly(
    size_t i,
    size_t j) const {
  FermiMonomial<Fermi1DLadderOp> mnI = OpSet[i];
  mnI.herm();
  FermiMonomial<Fermi1DLadderOp> mnJ = OpSet[j];
  mnI *= mnJ;
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > ans(mnI);
  ans.normalOrder();
  return ans;
}

//-------------------------------------------------------------Other Functions-----------

void printMatrixFermi1D(Fermi1DConsSet & constraints,
                        Fermi1DOpBasis & basis,
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
      FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > polyIJ =
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
  FermiTransMatToReIm(matrices, pairs);
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

void printSparseMatrixFermi(const Fermi1DConsSet & constraints,
                            const Fermi1DOpBasis & basis,
                            const std::string fileName,
                            const std::vector<std::complex<double> > ham,
                            const std::vector<std::pair<size_t, size_t> > & pairs,
                            bool isInf) {
  size_t matrixNum = basis.getLength();
  size_t matrixSize = constraints.getLength();
  vector<vector<vector<complex<double> > > > matrices(
      matrixNum,
      vector<vector<complex<double> > >(matrixSize,
                                        vector<complex<double> >(matrixSize)));
  for (size_t i = 0; i < matrixSize; i++) {
    for (size_t j = 0; j < matrixSize; j++) {
      FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > polyIJ =
          constraints.getIJPoly(i, j);
      //if ((polyIJ.getSize() % 2) != 0) {
      //  continue;
      //}
      //vector<complex<double> > entryIJ = basis.projPoly(polyIJ);
      try {
        vector<complex<double> > entryIJ = basis.projPolyInf(polyIJ);
        for (size_t k = 0; k < matrixNum; k++) {
          matrices[k][i][j] = entryIJ[k];
        }
      }
      catch (std::exception & e) {
      }
    }
  }
  FermiTransMatToReIm(matrices, pairs);
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
    ////////////////////////////////////
    else if (i == 2) {
      inputFile << -ham[i].real();
    }
    ///////////////////////////////////
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

void FermiTransMatToReIm(vector<vector<vector<complex<double> > > > & matrices,
                         const vector<pair<size_t, size_t> > & pairs) {
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

#endif  //QM_FERMI_CONSTRAINTS_NONTEM_CPP
