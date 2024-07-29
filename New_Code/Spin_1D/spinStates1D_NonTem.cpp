/*
  Jiazheng Sun
  Updated: Jul 28, 2024
  
  Class Implementations:
  SpinHalfBaseState1D
  SpinHalfState1D
  SpinHalfBasis1D
*/

#ifndef QM_SPINSTATES_1D_NONTEM_CPP
#define QM_SPINSTATES_1D_NONTEM_CPP

#include <cstddef>

#include "./spinStates1D.hpp"

using std::vector;

//----------------------------------------------------------SpinHalfBaseState1D----------

unsigned long SpinHalfBaseState1D::toDecimal() const {
  unsigned long decimalValue = 0;
  unsigned long base = 1;
  for (auto it = Nums.rbegin(); it != Nums.rend(); ++it) {
    if (*it) {
      decimalValue += base;
    }
    base <<= 1;  // Efficiently multiply base by 2
  }
  return decimalValue;
}

bool SpinHalfBaseState1D::operator<(const SpinHalfBaseState1D & rhs) const {
  return this->toDecimal() < rhs.toDecimal();
}

//-------------------------------------------------------------SpinHalfState1D-----------

SpinHalfState1D & SpinHalfState1D::operator=(const SpinHalfState1D & rhs) {
  if (this != &rhs) {
    Terms = rhs.Terms;
  }
  return *this;
}

//-------------------------------------------------------------SpinHalfBasis1D-----------

size_t VectorBoolHash::operator()(const SpinHalfBaseState1D & baseState) const {
  size_t hash = 0;
  std::vector<bool> vec = baseState.getAllNums();
  for (bool b : vec) {
    hash = (hash << 1) ^ std::hash<bool>()(b);
  }
  return hash;
}

void SpinHalfBasis1D::init() {
  size_t total = std::pow(2, SitesNum);
  for (size_t i = 0; i < total; i++) {
    std::bitset<32> bits(i);
    vector<bool> newState(SitesNum);
    for (size_t j = 0; j < SitesNum; j++) {
      newState[SitesNum - j - 1] = bits[j];
    }
    SpinHalfBaseState1D toAdd(newState);
    States.push_back(toAdd);
    lookupTable[toAdd] = States.size() - 1;
  }
}

void SpinHalfBasis1D::init(int SzTotal) {
  size_t total = std::pow(2, SitesNum);
  int requiredSpinUpCount = (SitesNum + SzTotal) / 2;
  for (size_t i = 0; i < total; i++) {
    int popCount = __builtin_popcount(i);
    if (popCount != requiredSpinUpCount) {
      continue;
    }
    vector<bool> newState(SitesNum);
    for (size_t j = 0; j < SitesNum; j++) {
      newState[SitesNum - j - 1] = ((i & (1 << j)) != 0);
    }
    SpinHalfBaseState1D toAdd(newState);
    States.push_back(toAdd);
    lookupTable[toAdd] = States.size() - 1;
  }
}

void SpinHalfBasis1D::initRefSym(int SzTotal) {
  size_t total = std::pow(2, SitesNum);
  int requiredSpinUpCount = (SitesNum + SzTotal) / 2;
  for (size_t i = 0; i < total; i++) {
    int popCount = __builtin_popcount(i);
    if (popCount != requiredSpinUpCount) {
      continue;
    }
    if (!isReverseEqual(i, SitesNum)) {
      if (!isLessOrEqualToReverse(i, SitesNum)) {
        continue;
      }
    }
    vector<bool> newState(SitesNum);
    for (size_t j = 0; j < SitesNum; j++) {
      newState[SitesNum - j - 1] = ((i & (1 << j)) != 0);
    }
    SpinHalfBaseState1D toAdd(newState);
    States.push_back(toAdd);
    lookupTable[toAdd] = States.size() - 1;
  }
}

std::string SpinHalfBasis1D::toString() {
  std::string ans = "SitesNum = ";
  ans += std::to_string(SitesNum);
  ans += "\nFull Basis:\n";
  for (vector<SpinHalfBaseState1D>::const_iterator it = States.begin();
       it != States.end();
       ++it) {
    ans += it->toString();
    ans += "\n";
  }
  return ans;
}

/*
int SpinHalfBasis1D::findBaseState(const SpinHalfBaseState1D & baseState) {
  return findIndex(States, baseState);
}
*/

size_t SpinHalfBasis1D::lookUpBaseState(const SpinHalfBaseState1D & baseState) {
  //return IndexTable[baseState];
  return lookupTable[baseState];
}

bool SpinHalfBasis1D::isReverseEqual(int num, int bitSize) {
  int reverseNum = 0;
  for (int i = 0; i < bitSize; ++i) {
    if (num & (1 << i)) {
      reverseNum |= (1 << (bitSize - 1 - i));
    }
  }
  return num == reverseNum;
}

bool SpinHalfBasis1D::isLessOrEqualToReverse(int num, int bitSize) {
  int reverseNum = 0;
  for (int i = 0; i < bitSize; ++i) {
    if (num & (1 << i)) {
      reverseNum |= (1 << (bitSize - 1 - i));
    }
  }
  return num <= reverseNum;
}

#endif  //QM_SPINSTATES_1D_NONTEM_CPP
