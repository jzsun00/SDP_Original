/*
  Jan 31, 2024  Jiazheng Sun

  Define ladder operators
*/

#ifndef ORI_SDP_GS_OPERATORS_HPP
#define ORI_SDP_GS_OPERATORS_HPP

#include <cstdlib>
#include <iostream>
#include <limits>

class LadderOp {
 protected:
  int index;
  bool creatorF;

 public:
  /*Constructs a ladder operator with specified index and creatorF,
   default constructor use INT_MIN and false.*/
  LadderOp() : index(std::numeric_limits<int>::min()), creatorF(false) {}
  LadderOp(int i_idx, bool i_creatorF) : index(i_idx), creatorF(i_creatorF) {}
  /*Get information of the ladder operator.*/
  int getIndex() const { return index; }
  bool getCreatorF() const { return creatorF; }
  /*Overload operators.*/
  bool operator==(LadderOp const & v) const {
    return creatorF == v.creatorF && index == v.index;
  }
  bool operator<(LadderOp const & v) const {
    if (creatorF != v.creatorF) {
      return creatorF < v.creatorF;
    }
    else {
      return index < v.index;
    }
  }
  bool operator>(LadderOp const & v) const { return !(*this < v || *this == v); }
  void herm() { creatorF ^= 1; }
};

#endif
