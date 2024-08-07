/*
  Jiazheng Sun
  Updated: Aug 7, 2024
*/

#ifndef QM_HARDCORE_CONSTRAINTS_HPP
#define QM_HARDCORE_CONSTRAINTS_HPP

#include "../Basics/constraints_Tem.hpp"
#include "./hardCoreOperators_Tem.hpp"
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
  virtual std::string toString() const;
};

//------------------------------------------------------------HardCore1DConsSet----------

class HardCore1DConsSet : public ConsSet<HardCoreMonomial<HardCore1DLadderOp>, int> {
 public:
  HardCore1DConsSet() : ConsSet<HardCoreMonomial<HardCore1DLadderOp>, int>() {
    HardCore1DLadderOp unit(true);
    OpSet.push_back(unit);
  }
  virtual std::string toString() const;
  virtual void addBaseSet(
      const ConsBaseSet<HardCoreMonomial<HardCore1DLadderOp>, int> & rhs);
  HardCorePolynomial<HardCoreMonomial<HardCore1DLadderOp> > getIJPoly(size_t i, size_t j);
};

//-------------------------------------------------------------Other Functions-----------

void printMatrixHardCore1D(HardCore1DConsSet & constraints,
                           HardCore1DOpBasis & basis,
                           std::string fileName,
                           std::vector<std::complex<double> > ham,
                           std::vector<std::pair<size_t, size_t> > & pairs);

void printMatrixXX1D(size_t max, std::string fileName);

void printSparseMatrixHardCore1D(HardCore1DConsSet & constraints,
                                 HardCore1DOpBasis & basis,
                                 std::string fileName,
                                 std::vector<std::complex<double> > ham,
                                 std::vector<std::pair<size_t, size_t> > & pairs);

void printSparseMatrixXX1D(size_t max, std::string fileName);

void transMatToReIm(
    std::vector<std::vector<std::vector<std::complex<double> > > > & matrices,
    std::vector<std::pair<size_t, size_t> > & pairs);

#endif  //QM_HARDCORE_CONSTRAINTS_HPP
