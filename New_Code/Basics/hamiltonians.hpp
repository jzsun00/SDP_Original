/*
  Jiazheng Sun
  Updated: Jul 26, 2024

  Class:
  SparseHamiltonian<PolyType, BaseStateType>
  SparseRealHamiltonian<PolyType, BaseStateType>
  FullHamiltonian<PolyType, BaseStateType>

  Define Hamiltonian matrices.
  SparseRealHamiltonian contains the lower triangular part of the matrix,
  since Hamiltonian matrices are Hermitian, real Hamiltonian matrices are symmetric.
  SparseHamiltonian and SparseRealHamiltonian can be used for ARPACK++,
  FullHamiltonian is mainly for testing and verification.
*/

#ifndef QM_HAMILTONIANS_HPP
#define QM_HAMILTONIANS_HPP

#include <cstddef>

#include "./operators_Tem.hpp"
#include "./states_Tem.hpp"

//------------------------------------SparseHamiltonian<PolyType, BaseStateType>---------

template<typename PolyType, typename BaseStateType>
class SparseHamiltonian {
 protected:
  PolyType poly;
  size_t dim;
  int nnz;                                   // Number of non-zero elements
  std::vector<int> irow;                     // Indices of row
  std::vector<int> pcol;                     // Positions for splitting
  std::vector<std::complex<double> > nzVal;  // Non-zero values

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
  std::vector<int> getIrow() const { return irow; }
  int getIrow(size_t i) const { return irow[i]; }
  std::vector<int> getPcol() const { return pcol; }
  int getPcol(size_t i) const { return pcol[i]; }
  std::vector<std::complex<double> > getNzVal() const { return nzVal; }
  std::complex<double> getNzVal(size_t i) const { return nzVal[i]; }
  std::string toString();
  /*Use the specified basis to create matrix.*/
  virtual void createMatrix(Basis<BaseStateType> & basis);
};

//----------------------------------SparseRealHamiltonian<PolyType, BaseStateType>-------

template<typename PolyType, typename BaseStateType>
class SparseRealHamiltonian {
 protected:
  PolyType poly;              // Hamiltonian operator
  size_t dim;                 // Dimension of the matrix
  size_t nnz;                 // Number of non-zero elements
  std::vector<int> irow;      // Indices of row
  std::vector<int> pcol;      // Positions for splitting
  std::vector<double> nzVal;  // Non-zero values

 public:
  /*Construct a sparse Hamiltonian with specified polynomial and dimension,
    default constructor use empty Hamiltonian.*/
  SparseRealHamiltonian() : poly(), dim(0), nnz(0), irow(), pcol(), nzVal() {}
  SparseRealHamiltonian(PolyType poly, size_t dim) :
      poly(poly), dim(dim), nnz(0), irow(), pcol(), nzVal() {}
  virtual ~SparseRealHamiltonian() {}
  /*Get information of the sparse Hamiltonian.*/
  size_t getDimension() const { return dim; }
  size_t getNumNonZero() const { return nnz; }
  std::vector<int> getIrow() const { return irow; }
  int getIrow(size_t i) const { return irow[i]; }
  std::vector<int> getPcol() const { return pcol; }
  int getPcol(size_t i) const { return pcol[i]; }
  std::vector<double> getNzVal() const { return nzVal; }
  std::complex<double> getNzVal(size_t i) const { return nzVal[i]; }
  std::string toString();
  /*Use the specified basis to create matrix.*/
  virtual void createMatrix(Basis<BaseStateType> & basis);
};

//-------------------------------------FullHamiltonian<PolyType, BaseStateType>----------

template<typename PolyType, typename BaseStateType>
class FullHamiltonian {
 protected:
  PolyType poly;
  size_t dim;
  std::vector<std::vector<std::complex<double> > > matrix;
  unsigned long nnz;

 public:
  /*Construct a full Hamiltonian, corresponding matrix is dim x dim.*/
  FullHamiltonian() : poly(), dim(0), matrix(), nnz(0) {}
  FullHamiltonian(PolyType poly, size_t dim) : poly(poly), dim(dim), matrix(), nnz(0) {
    matrix.resize(dim, std::vector<std::complex<double> >(dim, 0));
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

#endif  //QM_HAMILTONIANS_HPP
