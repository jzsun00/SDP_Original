/*
  Jiazheng Sun
  Updated: Aug 8, 2024
  
  Class Implementations:
  Fermi1DConsBaseSet
  Fermi1DConsSet
  
  Function Implementations:
*/

#ifndef QM_FERMI_CONSTRAINTS_NONTEM_CPP
#define QM_FERMI_CONSTRAINTS_NONTEM_CPP

#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>

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

//-------------------------------------------------------------Other Functions--------

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

void FermiPrintSparseSDPData(const Fermi1DConsSet & constraints,
                             const Fermi1DOpBasis & basis,
                             const std::string fileName,
                             const std::vector<std::complex<double> > & ham,
                             const std::vector<std::pair<size_t, size_t> > & pairs,
                             bool isInf) {
  size_t matrixNum = basis.getLength();
  size_t matrixSize = constraints.getLength();
  vector<ComplexCOOMatrix> COOMatrices(matrixNum,
                                       ComplexCOOMatrix(matrixSize, matrixSize));
  std::cout << "\nStart Matrix Construction" << std::endl;
  for (size_t i = 0; i < matrixSize; i++) {  //Fill the matrices
    for (size_t j = 0; j < matrixSize; j++) {
      FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > polyIJ =
          constraints.getIJPoly(i, j);
      vector<complex<double> > entryIJ(matrixNum);
      vector<size_t> validIdx;  //Only consider validIdx when adding to matrices
      try {                     //Project polyIJ using the basis
        if (isInf) {
          //entryIJ = basis.projPolyInf(polyIJ);
          validIdx = basis.projPolyInf(entryIJ, polyIJ);
        }
        else {
          validIdx = basis.projPolyFinite(entryIJ, polyIJ);
        }
        //std::cout << "polyIJ = " << polyIJ.toString() << std::endl;
        for (size_t k = 0; k < validIdx.size(); k++) {
          COOMatrices[validIdx[k]].addData(i, j, entryIJ[validIdx[k]]);
        }
      }
      catch (std::exception & e) {
        //std::cout << e.what() << std::endl;
      }
    }
  }
  FermiTransSparseMatToReIm(COOMatrices, pairs);
  COOMatrices[0] *= complex<double>(-1.0, 0);
  std::cout << "\nMatrix construction completed" << std::endl
            << "Start writing file" << std::endl;
  FermiPrintFileSparseSDPData(COOMatrices, ham, fileName);
}

void FermiPrintFileSparseSDPData(const vector<ComplexCOOMatrix> & COOMatrices,
                                 const std::vector<std::complex<double> > ham,
                                 const std::string & fileName) {
  std::ofstream inputFile(fileName);
  if (!inputFile.is_open()) {
    throw std::invalid_argument("ERROR: Failed to open file for writing!\n");
  }
  size_t matrixNum = COOMatrices.size();
  size_t matrixSize = COOMatrices[0].getNrows();
  inputFile << "\"XXZ Test: mDim = " << matrixNum - 1 << ", nBLOCK = 1, {"
            << matrixSize * 2 << "}\"" << std::endl;
  inputFile << matrixNum - 1 << "  =  mDIM" << std::endl;
  inputFile << "1  =  nBLOCK" << std::endl;
  inputFile << matrixSize * 2 << "  = bLOCKsTRUCT" << std::endl;
  inputFile << "{";
  for (size_t i = 1; i < ham.size(); i++) {  //Print cost function
    inputFile << ham[i].real();
    if (i < ham.size() - 1) {
      inputFile << ", ";
    }
  }
  inputFile << " }" << std::endl;
  for (size_t num = 0; num < matrixNum; num++) {  //Print constraint matrices
    size_t nnzNum = COOMatrices[num].getNnz();
    vector<size_t> rows = COOMatrices[num].getRows();
    vector<size_t> cols = COOMatrices[num].getCols();
    vector<complex<double> > values = COOMatrices[num].getAllData();
    for (size_t i = 0; i < nnzNum; ++i) {
      size_t row = rows[i];
      size_t col = cols[i];
      complex<double> value = values[i];
      // Real part
      if (row <= col && std::abs(value.real()) > ERROR) {
        inputFile << num << " "
                  << "1 " << row + 1 << " " << col + 1 << " " << value.real()
                  << std::endl;

        inputFile << num << " "
                  << "1 " << row + 1 + matrixSize << " " << col + 1 + matrixSize << " "
                  << value.real() << std::endl;
      }
      // Negative imaginary part (top right block)
      if (std::abs(value.imag()) > ERROR) {
        inputFile << num << " "
                  << "1 " << row + 1 << " " << col + 1 + matrixSize << " "
                  << -value.imag() << std::endl;
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

void FermiTransSparseMatToReIm(std::vector<ComplexCOOMatrix> & matrices,
                               const std::vector<std::pair<size_t, size_t> > & pairs) {
  size_t len = pairs.size();
  for (size_t n = 0; n < len; n++) {
    ComplexCOOMatrix mat1cp(matrices[pairs[n].first]);
    matrices[pairs[n].first] += matrices[pairs[n].second];
    matrices[pairs[n].second] *= complex<double>(-1.0, 0);
    matrices[pairs[n].second] += mat1cp;
    matrices[pairs[n].second] *= complex<double>(0, 1.0);
  }
}

#endif  //QM_FERMI_CONSTRAINTS_NONTEM_CPP
