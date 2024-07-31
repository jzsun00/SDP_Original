/*
  Jiazheng Sun
  Updated: Jul 31, 2024
*/

#ifndef QM_FERMI_SUBSPACES_HPP
#define QM_FERMI_SUBSPACES_HPP

#include <set>

#include "../Basics/subspaces_Tem.hpp"
#include "./fermiOperators_Tem.hpp"

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
  bool isNew(HardCoreMonomial<HardCore1DLadderOp> const & mn);
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
  vector<complex<double> > projPolyInf(
      HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > poly);
};

//--------------------------------------------------------------Other Functions----------

vector<pair<size_t, size_t> > findHermPairs(HardCore1DOpBasis & basis);

std::string printHermPairs(vector<pair<size_t, size_t> > & pairs);

void transVecToReIm(vector<complex<double> > & vec,
                    vector<pair<size_t, size_t> > & pairs);

#endif  //QM_FERMI_SUBSPACES_HPP
