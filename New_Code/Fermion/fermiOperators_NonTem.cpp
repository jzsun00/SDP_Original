/*
  Jiazheng Sun
  Updated: Aug 6, 2024
  
  Class Implementations:
  Fermi1DLadderOp
*/

#ifndef QM_FERMI_OPERATORS_NONTEM_CPP
#define QM_FERMI_OPERATORS_NONTEM_CPP

#include "./fermiOperators.hpp"

//------------------------------------------------------------Fermi1DLadderOp------------

Fermi1DLadderOp & Fermi1DLadderOp::operator=(const Fermi1DLadderOp & rhs) {
  if (this != &rhs) {
    this->index = rhs.index;
    this->isUnit = rhs.isUnit;
    this->creatorF = rhs.creatorF;
  }
  return *this;
}

bool Fermi1DLadderOp::operator<(const FermiLadderOp<int> & rhs) const {
  if (this->isUnit && rhs.getIsUnit()) {  //If both are e, they are equal
    return false;
  }
  if (this->isUnit && !rhs.getIsUnit()) {  //If only *this is e, *this is smaller
    return true;
  }
  if (rhs.getIsUnit() && !this->isUnit) {  //If only rhs is e, rhs is smaller
    return false;
  }
  if (!this->creatorF && !rhs.getCreatorF()) {  //Both are annihilation
    return this->index < rhs.getIndex();
  }
  else if (this->creatorF && rhs.getCreatorF()) {  //Both are creation
    return this->index > rhs.getIndex();
  }
  else {  //Not same type
    return this->creatorF < rhs.getCreatorF();
  }
}

#endif  //QM_FERMI_OPERATORS_NONTEM_CPP
