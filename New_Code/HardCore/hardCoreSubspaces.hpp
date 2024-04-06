/*
  Jiazheng Sun
  Updated: Apr 6, 2024
*/

#ifndef ORI_SDP_GS_HARDCORESUBSPACES_HPP
#define ORI_SDP_GS_HARDCORESUBSPACES_HPP

#include "../Basics/subspaces.hpp"
#include "./hardCoreOperators.hpp"

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

#endif
