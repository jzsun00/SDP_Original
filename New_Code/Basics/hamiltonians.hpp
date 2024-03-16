/*
  Define Hamiltonian matrices.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_HPP
#define ORI_SDP_GS_HAMILTONIANS_HPP

#include "./operators.hpp"
#include "./states.hpp"

//-------------------------------------------------------------SparseHamiltonian---------

template<typename PolyType, typename BasisType>
class SparseHamiltonian {
 private:
  PolyType poly;
  size_t dim;
  unsigned long nnz;
  vector<unsigned long> irow;
  vector<unsigned long> pcol;
  vector<complex<double> > nzVal;

 public:
  SparseHamiltonian() : poly(), dim(0), nnz(0), irow(), pcol(), nzVal() {}
  SparseHamiltonian(PolyType poly, size_t dim) :
      poly(poly), dim(dim), nnz(0), irow(), pcol(), nzVal() {}
  ~SparseHamiltonian() {}
  size_t getDimension() const { return dim; }
  unsigned long getNumNonZero() const { return nnz; }
  std::string toString();
  void createMatrix(BasisType & basis);
};

//-------------------------------------------------------------FullHamiltonian-----------

template<typename PolyType, typename BasisType>
class FullHamiltonian {
 protected:
  PolyType poly;
  size_t dim;
  vector<vector<complex<double> > > matrix;

 public:
  FullHamiltonian() : poly(), dim(0), matrix() {}
  FullHamiltonian(PolyType poly, size_t dim) : poly(poly), dim(dim), matrix() {
    matrix.resize(dim, std::vector<complex<double> >(dim, 0));
  }
  ~FullHamiltonian() {}
  size_t getDimension() const { return dim; }
  //unsigned long getNumNonZero() const;
  std::string toString();
  void createMatrix(BasisType const & basis);
};

#include "./hamiltonians_Tem.cpp"

#endif
