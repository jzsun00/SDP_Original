/*
  Jiazheng Sun
  Updated: Aug 15, 2024
  
  Class:
  Fermi1DConsBaseSet
  Fermi1DConsSet
  
  Function:
*/

#ifndef QM_FERMI_CONSTRAINTS_HPP
#define QM_FERMI_CONSTRAINTS_HPP

#include "../Basics/constraints_Tem.hpp"
#include "../LinearAlgebra/printTools.hpp"
#include "../LinearAlgebra/sparseCOO_Tem.hpp"
#include "./fermiOperators_Tem.hpp"
#include "./fermiSubspaces.hpp"

//----------------------------------------------------------Fermi1DConsBaseSet--------

class Fermi1DConsBaseSet : public ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int> {
 public:
  Fermi1DConsBaseSet() : ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int>() {}
  Fermi1DConsBaseSet(int start, int end, size_t order) :
      ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int>(start, end, order) {}
  Fermi1DConsBaseSet(const Fermi1DConsBaseSet & rhs) :
      ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int>(rhs) {}
  virtual void init();
  virtual ~Fermi1DConsBaseSet() {}
  /*Get information of the constraint base set.*/
  virtual std::string toString() const;
  /*Overload operators.*/
  Fermi1DConsBaseSet & operator=(const Fermi1DConsBaseSet & rhs);
};

//------------------------------------------------------------Fermi1DConsSet----------

class Fermi1DConsSet : public ConsSet<FermiMonomial<Fermi1DLadderOp>, int> {
 public:
  Fermi1DConsSet();
  Fermi1DConsSet(const Fermi1DConsSet & rhs) :
      ConsSet<FermiMonomial<Fermi1DLadderOp>, int>(rhs) {}
  virtual ~Fermi1DConsSet() {}
  /*Get information of the Fermi 1D constraint set.*/
  virtual std::string toString() const;
  /*Overload operators.*/
  Fermi1DConsSet & operator=(const Fermi1DConsSet & rhs);
  virtual void addBaseSet(const ConsBaseSet<FermiMonomial<Fermi1DLadderOp>, int> & rhs);
  FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > getIJPoly(size_t i, size_t j) const;
};

//-------------------------------------------------------------Other Functions--------

void printMatrixFermi1D(Fermi1DConsSet & constraints,
                        Fermi1DOpBasis & basis,
                        std::string fileName,
                        std::vector<std::complex<double> > ham,
                        std::vector<std::pair<size_t, size_t> > & pairs);

void FermiPrintSparseSDPData(const Fermi1DConsSet & constraints,
                             const Fermi1DOpBasis & basis,
                             const std::string fileName,
                             const std::vector<std::complex<double> > & ham,
                             const std::vector<std::pair<size_t, size_t> > & pairs,
                             bool isInf);

void FermiPrintFileSparseSDPData(const vector<ComplexCOOMatrix> & COOMatrices,
                                 const std::vector<std::complex<double> > ham,
                                 const std::string & fileName);

void FermiTransMatToReIm(
    std::vector<std::vector<std::vector<std::complex<double> > > > & matrices,
    const std::vector<std::pair<size_t, size_t> > & pairs);

void FermiTransSparseMatToReIm(std::vector<ComplexCOOMatrix> & matrices,
                               const std::vector<std::pair<size_t, size_t> > & pairs);

#endif  //QM_FERMI_CONSTRAINTS_HPP
