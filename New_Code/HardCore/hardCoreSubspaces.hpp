/*
  Jiazheng Sun
  Updated: Apr 6, 2024
*/

#ifndef ORI_SDP_GS_HARDCORESUBSPACES_HPP
#define ORI_SDP_GS_HARDCORESUBSPACES_HPP

#include "../Basics/subspaces.hpp"
#include "./hardCoreOperators.hpp"

//-------------------------------------------------------------HardCore1DOpSubBasis------

class HardCore1DOpSubBasis
    : public OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int> {
 public:
  HardCore1DOpSubBasis() : OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int>() {}
  HardCore1DOpSubBasis(int start, int end, size_t order) :
      OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int>(start, end, order) {}
  HardCore1DOpSubBasis(OpSubBasis const & rhs) :
      OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int>(rhs) {}
  virtual void init();
  ~HardCore1DOpSubBasis() {}
  virtual std::string toString();
};

//---------------------------------------------------------------HardCore1DOpBasis-------

class HardCore1DOpBasis : public OpBasis<HardCoreMonomial<HardCore1DLadderOp>, int> {
 public:
  HardCore1DOpBasis() : OpBasis<HardCoreMonomial<HardCore1DLadderOp>, int>() {
    HardCore1DLadderOp unit(true);
    Basis.push_back(unit);
  }
  ~HardCore1DOpBasis() {}
  HardCoreMonomial<HardCore1DLadderOp> operator[](size_t num) { return Basis[num]; }
  virtual std::string toString();
  virtual void addSubspace(OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int> & rhs);
};

#endif
