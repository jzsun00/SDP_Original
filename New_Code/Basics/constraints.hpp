/*
  Jiazheng Sun
  Updated: Aug 6, 2024
  
  Class:
  ConsBaseSet<MonomialType, IndexType>
  ConsSet<MonomialType, IndexType>
  
  Function:
  void printMatrix(ConsSet & constraints, OpBasis & basis);
*/

#ifndef QM_CONSTRAINTS_HPP
#define QM_CONSTRAINTS_HPP

#include "./subspaces_Tem.hpp"

//---------------------------------------ConsBaseSet<MonomialType, IndexType>------------

template<typename MonomialType, typename IndexType>
class ConsBaseSet {
 protected:
  IndexType start;
  IndexType end;
  size_t order;
  std::vector<MonomialType> BaseOpSet;

 public:
  ConsBaseSet() : order(0), BaseOpSet() {}
  ConsBaseSet(IndexType start, IndexType end, size_t order) :
      start(start), end(end), order(order), BaseOpSet() {}
  ConsBaseSet(const ConsBaseSet & rhs) :
      start(rhs.start), end(rhs.end), order(rhs.order), BaseOpSet(rhs.BaseOpSet) {}
  virtual void init() = 0;
  virtual ~ConsBaseSet() {}
  /*Get information of the constraint base set.*/
  size_t getLength() const { return BaseOpSet.size(); }
  IndexType getStart() const { return start; }
  IndexType getEnd() const { return end; }
  size_t getOrder() const { return order; }
  std::vector<MonomialType> getFullBaseOpSet() const { return BaseOpSet; }
  virtual std::string toString() = 0;
  /*Overload operators.*/
  ConsBaseSet & operator=(const ConsBaseSet<MonomialType, IndexType> & rhs);
  MonomialType operator[](size_t n) const { return BaseOpSet[n]; }
};

//------------------------------------------ConsSet<MonomialType, IndexType>-------------

template<typename MonomialType, typename IndexType>
class ConsSet {
 protected:
  std::vector<MonomialType> OpSet;

 public:
  ConsSet() : OpSet() {}
  ConsSet(const ConsSet<MonomialType, IndexType> & rhs) : OpSet(rhs.OpSet) {}
  virtual ~ConsSet() {}
  /*Get information of the constraint set.*/
  size_t getLength() const { return OpSet.size(); }
  virtual std::string toString() = 0;
  /*Overload operators.*/
  ConsSet & operator=(const ConsSet<MonomialType, IndexType> & rhs);
  virtual void addBaseSet(ConsBaseSet<MonomialType, IndexType> & rhs) = 0;
  //virtual Polynomial<MonomialType> getIJPoly(size_t i, size_t j);
};

//-------------------------------------------------------Other Functions-----------------

template<typename MonomialType, typename IndexType>
void printMatrix(ConsSet<MonomialType, IndexType> & constraints,
                 OpBasis<MonomialType, IndexType> & basis);

#endif  //QM_CONSTRAINTS_HPP
