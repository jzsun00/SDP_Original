/*
  Jiazheng Sun
  Updated: Jun 17, 2024

  Class:
  SpinHalfBaseState1D
  SpinHalfState1D
  SpinHalfBasis1D

  Define base states, states and basis for 1D spin-1/2 systems.
  Define the full basis for a 1D spin chain system.
*/

#ifndef QM_SPINSTATES1D_HPP
#define QM_SPINSTATES1D_HPP

#include <bitset>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <string>

#include "../Basics/states.hpp"
#include "../Basics/states_Tem.cpp"

//----------------------------------------------------------SpinHalfBaseState1D----------

class SpinHalfBaseState1D : public SpinBaseState<bool> {
 public:
  /*Construct a base state for 1D spin-1/2 system.
    Bool value 'true' for spin-up, 'false' for spin-down.
    Added constructor: if passing in size_t N, construct state with all spin-down.*/
  SpinHalfBaseState1D() : SpinBaseState() {}
  SpinHalfBaseState1D(size_t N) : SpinBaseState() { Nums = vector<bool>(N, false); }
  SpinHalfBaseState1D(vector<bool> & input) : SpinBaseState(input) {}
  SpinHalfBaseState1D(const SpinBaseState & rhs) : SpinBaseState(rhs) {}
  virtual ~SpinHalfBaseState1D() {}
  /*Get information of the spin base state.*/
  virtual std::string numToString(bool num) const { return std::to_string(num); }
  /*Overload operators.*/
  bool operator<(const SpinHalfBaseState1D & rhs) const;
};

//-------------------------------------------------------------SpinHalfState1D-----------

class SpinHalfState1D : public State<SpinHalfBaseState1D> {
 public:
  /*Construct a general quantum state for 1D spin-1/2 systems.
    Constructors are identical to State*/
  SpinHalfState1D() : State() {}
  SpinHalfState1D(const SpinHalfBaseState1D & ffs) : State(ffs) {}
  SpinHalfState1D(complex<double> pref, const SpinHalfBaseState1D & ffs) :
      State(pref, ffs) {}
  SpinHalfState1D(const SpinHalfState1D & rhs) : State(rhs) {}
  virtual ~SpinHalfState1D() {}
  /*Overload operators.*/
  SpinHalfState1D & operator=(const SpinHalfState1D & rhs);
};

//-------------------------------------------------------------SpinHalfBasis1D-----------

class SpinHalfBasis1D : public Basis<SpinHalfBaseState1D> {
 protected:
  size_t SitesNum;
  map<SpinHalfBaseState1D, size_t> IndexTable;

 public:
  /*Construct the entire basis of spin-1/2 systems.*/
  SpinHalfBasis1D() : Basis(), SitesNum(0), IndexTable() {}
  SpinHalfBasis1D(size_t n) : Basis(), SitesNum(n), IndexTable() {}
  virtual ~SpinHalfBasis1D() {}
  /*Fill the basis.
    If pass in an int SzTotal (should divide by 2 for actual Sz),
    construct a basis in that SzTotal sector.*/
  virtual void init();
  void init(int SzTotal);
  /*Get information of the full basis.*/
  virtual std::string toString();
  int findBaseState(const SpinHalfBaseState1D & baseState);
  size_t lookUpBaseState(const SpinHalfBaseState1D & baseState);
};

#endif  //QM_SPINSTATES1D_HPP
