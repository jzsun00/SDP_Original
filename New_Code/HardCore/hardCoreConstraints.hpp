/*
  Jiazheng Sun
  Updated: Apr 16, 2024
*/

#ifndef ORI_SDP_GS_HARDCORECONSTRAINTS_HPP
#define ORI_SDP_GS_HARDCORECONSTRAINTS_HPP

#include "../Basics/constraints.hpp"
#include "./hardCoreOperators.hpp"
#include "./hardCoreSubspaces.hpp"

//----------------------------------------------------------HardCore1DConsBaseSet--------

class HardCore1DConsBaseSet
    : public ConsBaseSet<HardCoreMonomial<HardCore1DLadderOp>, int> {
 public:
  HardCore1DConsBaseSet() : ConsBaseSet<HardCoreMonomial<HardCore1DLadderOp>, int>() {}
  HardCore1DConsBaseSet(int start, int end, size_t order) :
      ConsBaseSet<HardCoreMonomial<HardCore1DLadderOp>, int>(start, end, order) {}
  HardCore1DConsBaseSet(ConsBaseSet<HardCoreMonomial<HardCore1DLadderOp>, int> & rhs) :
      ConsBaseSet<HardCoreMonomial<HardCore1DLadderOp>, int>(rhs) {}
  virtual void init();
  ~HardCore1DConsBaseSet() {}
  virtual std::string toString();
};

//------------------------------------------------------------HardCore1DConsSet----------

class HardCore1DConsSet : public ConsSet<HardCoreMonomial<HardCore1DLadderOp>, int> {
 public:
  virtual std::string toString();
  virtual void addBaseSet(ConsBaseSet<HardCoreMonomial<HardCore1DLadderOp>, int> & rhs);
};

//-------------------------------------------------------------Other Functions-----------

void printMatrixHardCore1D(HardCore1DConsSet & constraints, HardCore1DOpBasis & basis);

#endif  //ORI_SDP_GS_HARDCORECONSTRAINTS_HPP
