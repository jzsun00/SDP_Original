/*
  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_HPP
#define ORI_SDP_GS_HAMILTONIANS_HPP

#include "./fermiOperators.hpp"
#include "./fermiStates.hpp"

template<typename polyType, typename basisType>
class Hamiltonian {
 private:
  polyType poly;
  size_t dim;
  unsigned long nnz;
  vector<unsigned long> irow;
  vector<unsigned long> pcol;
  vector<complex<double> > nzVal;

 public:
  Hamiltonian() : poly(), dim(0), nnz(0), irow(), pcol(), nzVal() {}
  Hamiltonian(polyType poly, size_t dim) :
      poly(poly), dim(dim), nnz(0), irow(), pcol(), nzVal() {}
  ~Hamiltonian() {}
  size_t getDimension() const { return dim; }
  unsigned long getNumNonZero() const { return nnz; }
  std::string toString();
  void createMatrix(basisType & basis);
};

#include "./hamiltonians.cpp"

#endif
