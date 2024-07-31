/*
  Jiazheng Sun
  Updated: Jul 31, 2024
  
  Class Implementations:
  Fermi1DLadderOp
*/

#ifndef QM_FERMI_OPERATORS_NONTEM_CPP
#define QM_FERMI_OPERATORS_NONTEM_CPP

#include "fermiOperators.hpp"

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
  if (this->isUnit && rhs.getIsUnit()) {
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
