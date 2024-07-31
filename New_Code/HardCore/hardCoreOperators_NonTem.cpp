/*
  Jiazheng Sun
  Updated: Jul 31, 2024

  Implementations of methods in class:
  Fermi1DLadderOp, FermiMonomial, FermiPolynomial.
 */

#ifndef QM_HARDCORE_OPERATORS_NONTEM_CPP
#define QM_HARDCORE_OPERATORS_NONTEM_CPP

#include "./hardCoreOperators.hpp"

//---------------------------------------------------------------HardCore1DLadderOp------

bool HardCore1DLadderOp::operator<(LadderOp const & rhs) const {
  if (this->isUnit && rhs.getIsUnit()) {
    return false;
  }
  /*
  if (this->creatorF == rhs.getCreatorF()) {
    if (this->creatorF == false) {
      return this->index < rhs.getIndex();
    }
    else {
      return this->index > rhs.getIndex();
    }
  }
  else {
    return this->creatorF < rhs.getCreatorF();
  }
  */
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

#endif  //QM_HARDCORE_OPERATORS_NONTEM_CPP
