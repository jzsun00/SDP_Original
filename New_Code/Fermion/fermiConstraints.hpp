/*
  Jiazheng Sun
  Updated: Aug 2, 2024
*/

#ifndef QM_FERMI_CONSTRAINTS_HPP
#define QM_FERMI_CONSTRAINTS_HPP

#include "../Basics/constraints_Tem.hpp"
#include "./fermiOperators_Tem.hpp"
#include "./fermiSubspaces.hpp"

//----------------------------------------------------------Fermi1DConsBaseSet--------

class Fermi1DConsBaseSet : public ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int> {
 public:
  Fermi1DConsBaseSet() : ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int>() {}
  Fermi1DConsBaseSet(int start, int end, size_t order) :
      ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int>(start, end, order) {}
  Fermi1DConsBaseSet(ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int> & rhs) :
      ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int>(rhs) {}
  virtual void init();
  virtual ~Fermi1DConsBaseSet() {}
  virtual std::string toString();
};

//------------------------------------------------------------Fermi1DConsSet----------

class Fermi1DConsSet : public ConsSet<FermiMonomial<Fermi1DLadderOp>, int> {
 public:
  Fermi1DConsSet() : ConsSet<FermiMonomial<Fermi1DLadderOp>, int>() {
    Fermi1DLadderOp unit(true);
    OpSet.push_back(unit);
  }
  virtual ~Fermi1DConsSet() {}
  virtual std::string toString();
  virtual void addBaseSet(ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int> & rhs);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > getIJPoly(size_t i, size_t j);
};

//-------------------------------------------------------------Other Functions-----------

void printMatrixFermi1D(Fermi1DConsSet & constraints,
                        Fermi1DOpBasis & basis,
                        std::string fileName,
                        std::vector<std::complex<double> > ham,
                        std::vector<std::pair<size_t, size_t> > & pairs);

void printSparseMatrixFermi1D(Fermi1DConsSet & constraints,
                              Fermi1DOpBasis & basis,
                              std::string fileName,
                              std::vector<std::complex<double> > ham,
                              std::vector<std::pair<size_t, size_t> > & pairs);

void FermiTransMatToReIm(
    std::vector<std::vector<std::vector<std::complex<double> > > > & matrices,
    std::vector<std::pair<size_t, size_t> > & pairs);

#endif  //QM_FERMI_CONSTRAINTS_HPP
