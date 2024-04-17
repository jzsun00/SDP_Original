/*
  Jiazheng Sun
  Updated: Apr 16, 2024
*/

#ifndef ORI_SDP_GS_CONSTRAINTS_HPP
#define ORI_SDP_GS_CONSTRAINTS_HPP

#include "./operators.hpp"
#include "./subspaces.hpp"

//---------------------------------------------------------------ConsBaseSet-------------

template<typename MonomialType, typename IndexType>
class ConsBaseSet {
 protected:
  IndexType start;
  IndexType end;
  size_t order;
  vector<MonomialType> BaseOpSet;

 public:
  ConsBaseSet() : order(0), BaseOpSet() {}
  ConsBaseSet(IndexType start, IndexType end, size_t order) :
      start(start), end(end), order(order), BaseOpSet() {}
  ConsBaseSet(ConsBaseSet const * rhs) : BaseOpSet(rhs.BaseOpSet) {}
  virtual void init() = 0;
  ~ConsBaseSet() {}
  size_t getLength() const { return BaseOpSet.size(); }
  IndexType getStart() const { return start; }
  IndexType getEnd() const { return end; }
  size_t getOrder() const { return order; }
  vector<MonomialType> getBaseOpSet() const { return BaseOpSet; }
  MonomialType operator[](size_t n) const { return BaseOpSet[n]; }
  virtual std::string toString() = 0;
};

//-----------------------------------------------------------------ConsSet---------------

template<typename MonomialType, typename IndexType>
class ConsSet {
 protected:
  vector<MonomialType> OpSet;

 public:
  ConsSet() : OpSet() {}
  ~ConsSet() {}
  size_t getLength() const { return OpSet.size(); }
  virtual std::string toString() = 0;
  virtual void addBaseSet(ConsBaseSet<MonomialType, IndexType> & rhs) = 0;
  //virtual Polynomial<MonomialType> getIJPoly(size_t i, size_t j);
};

//-------------------------------------------------------------Other Functions-----------

template<typename MonomialType, typename IndexType>
void printMatrix(ConsSet<MonomialType, IndexType> & constraints,
                 OpBasis<MonomialType, IndexType> & basis);

#include "./constraints_Tem.cpp"

#endif  //ORI_SDP_GS_CONSTRAINTS_HPP
