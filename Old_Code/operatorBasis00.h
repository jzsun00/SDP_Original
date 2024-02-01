//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 29.06.2011 Thomas Barthel
//% SDP groundstate kit
//% Generating the operator basis for a certain choice of subsystems.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// key words: TODO, QUESTION, TEST, CHECK, REMARK, CAUTION
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#ifndef SDP_GS_KIT_OPERATORBASIS00_H
#define SDP_GS_KIT_OPERATORBASIS00_H

#include "algebra/operator00.h"
#include "lattice/Subsystems00.h"



//__________________________________________________________________________________________________
//###################################################################################### FockState #
class FockState
 {
 public:
 typedef LinearAlgebra::Vector<unsigned long> container_type;
 
 FockState(size_t i_N, size_t i_NoModes, bool i_PauliF) : pauliF(i_PauliF), Nmodes(i_NoModes), validF(1)
  {
  if( ( pauliF && i_NoModes < i_N ) || ( !pauliF && !i_NoModes ) )  validF = 0;
  if( pauliF )  idx = LinearAlgebra::range(0,i_N);
  else          idx = container_type(i_N,0);
  }
 
 size_t NoParticles() const  { return idx.size(); }
 size_t NoModes()     const  { return Nmodes; }
 bool   PauliF()      const  { return pauliF; } 
 bool   is_valid()    const  { return validF; }
 container_type const& operator()() const  { return idx; }
 
 inline bool inc();
 
 protected:
 bool           pauliF;
 size_t         Nmodes;
 bool           validF;
 container_type idx;
 };


//--------------------------------------------------------------------------------[ inc ]-
inline bool FockState::inc()
 {
 for(size_t n = 0; n < NoParticles(); n++)
  {
  if( ( ( n+1 < NoParticles() && idx[n]+1+pauliF <= idx[n+1] ) || n+1 == NoParticles() ) && idx[n]+1 < Nmodes )
   {
   idx[n]++;
   return true;
   }
  else
   {
   if( pauliF )  idx[LinearAlgebra::range(0,n+1)] = LinearAlgebra::range(0,n+1);
   else          idx[LinearAlgebra::range(0,n+1)] = container_type(n+1,0);
   }
  }
 validF = 0;
 return false;
 }



//##################################################################################################

//--------------------------------------------------------------------------[ factorial ]-
double log_factorial(unsigned n);

//---------------------------------------------------------------------------[ binomial ]-
double    log_binomial(unsigned D, unsigned n);
inline size_t binomial(unsigned D, unsigned n)  { return (size_t)floor( exp(log_binomial(D,n)) + 0.5 ); }


//-----------------------------------------------------------------------[ NoFockStates ]-
size_t NoFockStates(size_t N, size_t NoModes, bool PauliF);

inline size_t NoFockStates(FockState const& v)
 { return NoFockStates( v.NoParticles(), v.NoModes(), v.PauliF() ); }



//##################################################################################################

//-----------------------------------------------------------------[ operatorBasis_full ]-
// On exit, NoElements[K] gives the number of basis elements that are (k<=K)-point monomials,
//  giving, at the same time, the position of the first (K+1)-point term in the returned basis.
Monomials operatorBasis_full(Subsystems const& Omega, int algebraF, sdp_gs::Indices& NoElements);

inline Monomials operatorBasis_full(Subsystems const& Omega, int algebraF)
 {
 sdp_gs::Indices NoElements;
 return operatorBasis_full( Omega, algebraF, NoElements );
 }



//##################################################################################################

#endif // SDP_GS_KIT_OPERATORBASIS00_H
