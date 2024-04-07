/*
  Jiazheng Sun
  Updated: Apr 6, 2024
*/

#ifndef ORI_SDP_GS_SUBSPACES_HPP
#define ORI_SDP_GS_SUBSPACES_HPP

#include "./operators.hpp"

//---------------------------------------------------------------OpSubBasis--------------

template<typename MonomialType, typename IndexType>
class OpSubBasis {
 protected:
  IndexType start;
  IndexType end;
  size_t order;
  vector<MonomialType> Basis;

 public:
  /*Construct a basis with MonomialType.
    Default constructor use random size and empty vector.*/
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
  vector<MonomialType> getBasis() const { return Basis; }
  MonomialType operator[](size_t n) const { return Basis[n]; }
  virtual std::string toString() = 0;

  /*Projection tools.*/
  vector<complex<double> > projPoly(Polynomial<MonomialType> poly);
};

//---------------------------------------------------------------OpBasis-----------------

template<typename MonomialType, typename IndexType>
class OpBasis {
 protected:
  vector<MonomialType> Basis;

 public:
  OpBasis() : Basis() {}
  ~OpBasis() {}
  size_t getLength() const { return Basis.size(); }
  virtual std::string toString() = 0;
  virtual void addSubspace(OpSubBasis<MonomialType, IndexType> & rhs) = 0;

  /*Projection tools.*/
  vector<complex<double> > projPoly(Polynomial<MonomialType> poly);
};

#include "./subspaces_Tem.hpp"

#endif  //ORI_SDP_GS_SUBSPACES_HPP
