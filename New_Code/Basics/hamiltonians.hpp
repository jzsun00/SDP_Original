/*
  Jiazheng Sun
  Updated: Mar 17, 2024

  Class:
  Sparsehamiltonian, FullHamiltonian.

  Define Hamiltonian matrices.
  SparseHamiltonian can be used for ARPACK++,
  FullHamiltonian is mainly for testing and verification.
*/

#ifndef ORI_SDP_GS_HAMILTONIANS_HPP
#define ORI_SDP_GS_HAMILTONIANS_HPP

#include "./operators.hpp"
#include "./states.hpp"

//-------------------------------------------------------------SparseHamiltonian---------

template<typename PolyType, typename BaseStateType>
class SparseHamiltonian {
 private:
  PolyType poly;
  size_t dim;
  int nnz;                         // Number of non-zero elements
  vector<int> irow;                // Indices of row
  vector<int> pcol;                // Positions for splitting
  vector<complex<double> > nzVal;  // Non-zero values

 public:
  /*Construct a sparse Hamiltonian with specified polynomial and dimension,
    default constructor use empty Hamiltonian.*/
  SparseHamiltonian() : poly(), dim(0), nnz(0), irow(), pcol(), nzVal() {}
  SparseHamiltonian(PolyType poly, size_t dim) :
      poly(poly), dim(dim), nnz(0), irow(), pcol(), nzVal() {}
  ~SparseHamiltonian() {}
  /*Get information of the sparse Hamiltonian.*/
  size_t getDimension() const { return dim; }
  int getNumNonZero() const { return nnz; }
  vector<int> getIrow() const { return irow; }
  int getIrow(size_t i) const { return irow[i]; }
  vector<int> getPcol() const { return pcol; }
  int getPcol(size_t i) const { return pcol[i]; }
  vector<complex<double> > getNzVal() const { return nzVal; }
  complex<double> getNzVal(size_t i) const { return getNzVal[i]; }
  std::string toString();
  /*Use the specified basis to create matrix.*/
  void createMatrix(Basis<BaseStateType> & basis);
};

//-------------------------------------------------------------FullHamiltonian-----------

template<typename PolyType, typename BaseStateType>
class FullHamiltonian {
 protected:
  PolyType poly;
  size_t dim;
  vector<vector<complex<double> > > matrix;
  unsigned long nnz;

 public:
  /*Construct a full Hamiltonian, corresponding matrix is dim x dim.*/
  FullHamiltonian() : poly(), dim(0), matrix(), nnz(0) {}
  FullHamiltonian(PolyType poly, size_t dim) : poly(poly), dim(dim), matrix(), nnz(0) {
    matrix.resize(dim, std::vector<complex<double> >(dim, 0));
  }
  ~FullHamiltonian() {}
  /*Get information of the full Hamiltonian.*/
  PolyType getPoly() const { return poly; }
  size_t getDimension() const { return dim; }
  unsigned long getNumNonZero() const { return nnz; };
  std::string toString();
  /*Create the corresponding matrix in the given basis.*/
  void createMatrix(Basis<BaseStateType> const & basis);
};

#include "./hamiltonians_Tem.cpp"

#endif  //ORI_SDP_GS_HAMILTONIANS_HPP
