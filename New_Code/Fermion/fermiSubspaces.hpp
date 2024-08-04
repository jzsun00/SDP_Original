/*
  Jiazheng Sun
  Updated: Aug 3, 2024
  
  Class:
  Fermi1DOpSubBasis
  Fermi1DOpBasis
  
  Function:
  vector<pair<size_t, size_t> > findHermPairs(const Fermi1DOpBasis & basis)
  string printHermPairs(const vector<pair<size_t, size_t> > & pairs)
  void transVecToReIm(vector<complex<double> > & vec, vector<pair<size_t, size_t> > & pairs)
*/

#ifndef QM_FERMI_SUBSPACES_HPP
#define QM_FERMI_SUBSPACES_HPP

#include "../Basics/subspaces_Tem.hpp"
#include "./fermiOperators_Tem.hpp"

//----------------------------------------------------------Fermi1DOpSubBasis------------

class Fermi1DOpSubBasis : public OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int> {
 public:
  Fermi1DOpSubBasis() : OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int>() {}
  Fermi1DOpSubBasis(int start, int end, size_t order) :
      OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int>(start, end, order) {}
  Fermi1DOpSubBasis(const OpSubBasis & rhs) :
      OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int>(rhs) {}
  virtual void init();
  virtual ~Fermi1DOpSubBasis() {}
  /*Get information of the operator sub-basis.*/
  virtual std::string toString();
  bool isNew(const FermiMonomial<Fermi1DLadderOp> & toAdd) const;
};

//-----------------------------------------------------------Fermi1DOpBasis--------------

class Fermi1DOpBasis : public OpBasis<FermiMonomial<Fermi1DLadderOp>, int> {
 public:
  Fermi1DOpBasis() : OpBasis<FermiMonomial<Fermi1DLadderOp>, int>() {
    Fermi1DLadderOp unit(true);
    Basis.push_back(unit);
  }
  virtual ~Fermi1DOpBasis() {}
  FermiMonomial<Fermi1DLadderOp> operator[](size_t num) { return Basis[num]; }
  /*Get information of the operator basis.*/
  virtual std::string toString();
  virtual void addSubspace(const OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int> & rhs);
  std::vector<std::complex<double> > projPolyInf(
      FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > poly);
};

//-----------------------------------------------------------Other Functions-------------

std::vector<std::pair<size_t, size_t> > FermiFindHermPairs(Fermi1DOpBasis & basis);

std::string FermiPrintHermPairs(std::vector<std::pair<size_t, size_t> > & pairs);

void FermiTransVecToReIm(std::vector<std::complex<double> > & vec,
                         std::vector<std::pair<size_t, size_t> > & pairs);

#endif  //QM_FERMI_SUBSPACES_HPP
