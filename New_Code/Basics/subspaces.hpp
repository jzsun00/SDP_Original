/*
  Jiazheng Sun
  Updated: Jul 26, 2024
*/

#ifndef QM_SUBSPACES_HPP
#define QM_SUBSPACES_HPP

#include "./operators_Tem.hpp"

//---------------------------------------------------------------OpSubBasis--------------

template<typename MonomialType, typename IndexType>
class OpSubBasis {
 protected:
  IndexType start;
  IndexType end;
  size_t order;
  std::vector<MonomialType> Basis;

 public:
  /*Construct a basis with MonomialType.
    Default constructor use random size and empty std::vector.*/
  OpSubBasis() : order(0), Basis() {}
  OpSubBasis(IndexType start, IndexType end, size_t order) :
      start(start), end(end), order(order), Basis() {}
  OpSubBasis(OpSubBasis const & rhs) : Basis(rhs.Basis) {}
  virtual void init() = 0;
  ~OpSubBasis() {}
  /*Get information of the operator basis.*/
  size_t getLength() const { return Basis.size(); }
  IndexType getStart() const { return start; }
  IndexType getEnd() const { return end; }
  size_t getOrder() const { return order; }
  std::vector<MonomialType> getBasis() const { return Basis; }
  MonomialType operator[](size_t n) const { return Basis[n]; }
  virtual std::string toString() = 0;

  /*Projection tools.*/
  std::vector<std::complex<double> > projPoly(Polynomial<MonomialType> poly);
};

//---------------------------------------------------------------OpBasis-----------------

template<typename MonomialType, typename IndexType>
class OpBasis {
 protected:
  std::vector<MonomialType> Basis;

 public:
  OpBasis() : Basis() {}
  ~OpBasis() {}
  size_t getLength() const { return Basis.size(); }
  virtual std::string toString() = 0;
  virtual void addSubspace(OpSubBasis<MonomialType, IndexType> & rhs) = 0;

  /*Projection tools.*/
  std::vector<std::complex<double> > projPoly(Polynomial<MonomialType> poly);
};

#endif  //QM_SUBSPACES_HPP
