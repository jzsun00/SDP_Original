/*
  Jiazheng Sun
  Updated: Aug 7, 2024
  
  Class:
  Fermi1DOpSubBasis
  Fermi1DOpBasis
  
  Function:
  vector<pair<size_t, size_t> > findHermPairs(const Fermi1DOpBasis & basis)
  string printHermPairs(const vector<pair<size_t, size_t> > & pairs)
  void transVecToReIm(vector<complex<double> > & vec, vector<pair<size_t, size_t> > & pairs)
  
  Define basis and sub-basis for Fermi Operators.
*/

#ifndef QM_FERMI_SUBSPACES_HPP
#define QM_FERMI_SUBSPACES_HPP

#include "../Basics/subspaces_Tem.hpp"
#include "./fermiOperators.hpp"
#include "./fermiOperators_Tem.hpp"

//----------------------------------------------------------Fermi1DOpSubBasis------------

class Fermi1DOpSubBasis : public OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int> {
 public:
  Fermi1DOpSubBasis() : OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int>() {}
  Fermi1DOpSubBasis(int start, int end, size_t order) :
      OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int>(start, end, order) {}
  Fermi1DOpSubBasis(const Fermi1DOpSubBasis & rhs) :
      OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int>(rhs) {}
  virtual ~Fermi1DOpSubBasis() {}
  /*Overload operators.*/
  Fermi1DOpSubBasis & operator=(const Fermi1DOpSubBasis & rhs);
  virtual void init(bool isInf);
  /*Get information of the operator sub-basis.*/
  virtual std::string toString() const;
  bool isNew(const FermiMonomial<Fermi1DLadderOp> & toAdd) const;
};

//-----------------------------------------------------------Fermi1DOpBasis--------------

class Fermi1DOpBasis : public OpBasis<FermiMonomial<Fermi1DLadderOp>, int> {
 protected:
  std::map<FermiMonomial<Fermi1DLadderOp>, size_t> lookupTable;

 public:
  Fermi1DOpBasis();
  Fermi1DOpBasis(const Fermi1DOpBasis & rhs) :
      OpBasis<FermiMonomial<Fermi1DLadderOp>, int>(rhs) {}
  virtual ~Fermi1DOpBasis() {}
  /*Overload operators.*/
  Fermi1DOpBasis & operator=(const Fermi1DOpBasis & rhs);
  FermiMonomial<Fermi1DLadderOp> operator[](size_t num) const { return Basis[num]; }
  void buildTable();
  /*Get information of the operator basis.*/
  virtual std::string toString() const;
  virtual void addSubspace(const OpSubBasis<FermiMonomial<Fermi1DLadderOp>, int> & rhs);
  size_t findIndex(const FermiMonomial<Fermi1DLadderOp> & mn) const;
  std::vector<std::complex<double> > projPolyInf(
      const FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > & poly);
  std::vector<size_t> projPolyFinite(
      std::vector<std::complex<double> > & vec,
      const FermiPolynomial<FermiMonomial<Fermi1DLadderOp> > & poly);
};

//-----------------------------------------------------------Other Functions-------------

std::vector<std::pair<size_t, size_t> > FermiFindHermPairs(const Fermi1DOpBasis & basis);

std::string FermiPrintHermPairs(const std::vector<std::pair<size_t, size_t> > & pairs);

void FermiTransVecToReIm(std::vector<std::complex<double> > & vec,
                         const std::vector<std::pair<size_t, size_t> > & pairs);

#endif  //QM_FERMI_SUBSPACES_HPP
