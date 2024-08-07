/*
  Jiazheng Sun
  Updated: Aug 6, 2024
  
  Class:
  OpSubBasis<MonomialType, IndexType>
  OpBasis<MonomialType, IndexType>
*/

#ifndef QM_SUBSPACES_HPP
#define QM_SUBSPACES_HPP

#include "./operators_Tem.hpp"

//-----------------------------------------OpSubBasis<MonomialType, IndexType>-----------

template<typename MonomialType, typename IndexType>
class OpSubBasis {
 protected:
  IndexType start;  //Start site index
  IndexType end;    //End site index
  size_t order;     //Order of the Monomial
  std::vector<MonomialType> Basis;

 public:
  /*Construct a basis with MonomialType.
    Default constructor use random size and empty std::vector.*/
  OpSubBasis() : order(0), Basis() {}
  OpSubBasis(IndexType start, IndexType end, size_t order) :
      start(start), end(end), order(order), Basis() {}
  OpSubBasis(const OpSubBasis & rhs) :
      start(rhs.start), end(rhs.end), order(rhs.order), Basis(rhs.Basis) {}
  virtual void init() = 0;
  virtual ~OpSubBasis() {}
  /*Get information of the operator sub-basis.*/
  size_t getLength() const { return Basis.size(); }
  IndexType getStart() const { return start; }
  IndexType getEnd() const { return end; }
  size_t getOrder() const { return order; }
  std::vector<MonomialType> getFullBasis() const { return Basis; }
  virtual std::string toString() const = 0;
  /*Overload operators.*/
  OpSubBasis & operator=(const OpSubBasis<MonomialType, IndexType> & rhs);
  MonomialType operator[](size_t n) const { return Basis[n]; }
  /*Projection tools.*/
  virtual std::vector<std::complex<double> > projPoly(
      const Polynomial<MonomialType> & poly) const;
};

//-------------------------------------------OpBasis<MonomialType, IndexType>------------

template<typename MonomialType, typename IndexType>
class OpBasis {
 protected:
  std::vector<MonomialType> Basis;

 public:
  OpBasis() : Basis() {}
  OpBasis(const OpBasis<MonomialType, IndexType> & rhs) : Basis(rhs.Basis) {}
  virtual ~OpBasis() {}
  /*Get information of the operator basis.*/
  size_t getLength() const { return Basis.size(); }
  virtual std::string toString() const = 0;
  virtual void addSubspace(const OpSubBasis<MonomialType, IndexType> & rhs) = 0;
  /*Overload operators.*/
  OpBasis & operator=(const OpBasis<MonomialType, IndexType> & rhs);
  /*Projection tools.*/
  virtual std::vector<std::complex<double> > projPoly(
      const Polynomial<MonomialType> & poly) const;
};

#endif  //QM_SUBSPACES_HPP
