/*
  Jiazheng Sun
  Updated: Aug 15, 2024
  
  Class:
  SpinHalfBaseState1D
  SpinHalfState1D
  SpinHalfBasis1D
  
  Define base states, states and basis for 1D spin systems.
*/

#ifndef QM_SPIN_STATES_1D_HPP
#define QM_SPIN_STATES_1D_HPP

#include <bit>
#include <bitset>
#include <cmath>
#include <functional>
#include <numeric>
#include <unordered_map>

#include "../Basics/states_Tem.hpp"

//----------------------------------------------------------SpinHalfBaseState1D----------

class SpinHalfBaseState1D : public SpinBaseState<bool> {
 public:
  /*Construct a base state for 1D spin-1/2 system.
    Bool value 'true' for spin-up, 'false' for spin-down.
    Added constructor: if passing in size_t, construct state with all spin-down.*/
  SpinHalfBaseState1D() : SpinBaseState() {}
  SpinHalfBaseState1D(size_t N) : SpinBaseState(std::vector<bool>(N, false)) {}
  SpinHalfBaseState1D(const std::vector<bool> & input) : SpinBaseState(input) {}
  SpinHalfBaseState1D(const SpinHalfBaseState1D & rhs) : SpinBaseState(rhs) {}
  virtual ~SpinHalfBaseState1D() {}
  /*Get information of the 1D spin-1/2 base state.*/
  virtual std::string numToString(bool num) const { return std::to_string(num); }
  unsigned long toDecimal() const;  //Corresponding decimal value
  /*Overload operators.*/
  SpinHalfBaseState1D & operator=(const SpinHalfBaseState1D & rhs);
  bool operator<(const SpinHalfBaseState1D & rhs) const;  //Based on decimal value
};

//-------------------------------------------------------------SpinHalfState1D-----------

class SpinHalfState1D : public State<SpinHalfBaseState1D> {
 public:
  /*Construct a general quantum state for 1D spin-1/2 systems.
    Constructors are identical to State*/
  SpinHalfState1D() : State() {}
  SpinHalfState1D(const SpinHalfBaseState1D & rhs) : State(rhs) {}
  SpinHalfState1D(std::complex<double> pref, const SpinHalfBaseState1D & rhs) :
      State(pref, rhs) {}
  SpinHalfState1D(const SpinHalfState1D & rhs) : State(rhs) {}
  virtual ~SpinHalfState1D() {}
  /*Overload operators.*/
  SpinHalfState1D & operator=(const SpinHalfState1D & rhs);
  void clear();
};

//-------------------------------------------------------------SpinHalfBasis1D-----------

struct VectorBoolHash {
  size_t operator()(const SpinHalfBaseState1D & baseState) const;
};

class SpinHalfBasis1D : public Basis<SpinHalfBaseState1D> {
 protected:
  size_t SitesNum;  //Number of sites
  std::unordered_map<SpinHalfBaseState1D, size_t, VectorBoolHash>
      lookupTable;  //<BaseState, index> pair for searching

 public:
  /*Construct the entire basis of spin-1/2 systems.*/
  SpinHalfBasis1D() : Basis(), SitesNum(0), lookupTable() {}
  SpinHalfBasis1D(size_t n) : Basis(), SitesNum(n), lookupTable() {}
  SpinHalfBasis1D(const SpinHalfBasis1D & rhs);
  virtual ~SpinHalfBasis1D() {}
  /*Fill the basis.
    If pass in an int SzTotal (should divide by 2 for actual Sz),
    construct a basis in that SzTotal sector.*/
  virtual void init();
  void init(int SzTotal);
  void initRefSym(int SzTotal);
  /*Get information of the full basis.*/
  virtual std::string toString();
  int findBaseState(const SpinHalfBaseState1D & baseState);
  size_t lookUpBaseState(const SpinHalfBaseState1D & baseState);

 private:
  bool isReverseEqual(int num, int bitSize);
  bool isLessOrEqualToReverse(int num, int bitSize);
};

#endif  //QM_SPIN_STATES_1D_HPP
