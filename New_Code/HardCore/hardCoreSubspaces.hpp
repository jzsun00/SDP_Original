/*
  Jiazheng Sun
  Updated: Aug 6, 2024
*/

#ifndef QM_HARDCORE_SUBSPACES_HPP
#define QM_HARDCORE_SUBSPACES_HPP

#include <complex>
#include <set>

#include "../Basics/subspaces_Tem.hpp"
#include "./hardCoreOperators_Tem.hpp"

//-------------------------------------------------------------HardCore1DOpSubBasis------

class HardCore1DOpSubBasis
    : public OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int> {
 public:
  HardCore1DOpSubBasis() : OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int>() {}
  HardCore1DOpSubBasis(int start, int end, size_t order) :
      OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int>(start, end, order) {}
  HardCore1DOpSubBasis(HardCore1DOpSubBasis const & rhs) :
      OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int>(rhs) {}
  virtual void init(bool isInf);
  virtual ~HardCore1DOpSubBasis() {}
  virtual std::string toString() const;
  bool isNew(HardCoreMonomial<HardCore1DLadderOp> const & mn);
};

//---------------------------------------------------------------HardCore1DOpBasis-------

class HardCore1DOpBasis : public OpBasis<HardCoreMonomial<HardCore1DLadderOp>, int> {
 public:
  HardCore1DOpBasis() : OpBasis<HardCoreMonomial<HardCore1DLadderOp>, int>() {
    HardCore1DLadderOp unit(true);
    Basis.push_back(unit);
  }
  virtual ~HardCore1DOpBasis() {}
  HardCoreMonomial<HardCore1DLadderOp> operator[](size_t num) { return Basis[num]; }
  virtual std::string toString() const;
  virtual void addSubspace(
      const OpSubBasis<HardCoreMonomial<HardCore1DLadderOp>, int> & rhs);
  std::vector<std::complex<double> > projPolyInf(
      HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly);
};

//--------------------------------------------------------------Other Functions----------

std::vector<std::pair<size_t, size_t> > findHermPairs(HardCore1DOpBasis & basis);

std::string printHermPairs(std::vector<std::pair<size_t, size_t> > & pairs);

void transVecToReIm(std::vector<std::complex<double> > & vec,
                    std::vector<std::pair<size_t, size_t> > & pairs);

#endif  //ORI_SDP_GS_HARDCORESUBSPACES_HPP
