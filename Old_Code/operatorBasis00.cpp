#include "operatorBasis00.h"


//##################################################################################################

//--------------------------------------------------------------------------[ factorial ]-
double log_factorial(unsigned n)
 {
 double ret=0;
 for(unsigned m=2; m<=n; m++)  ret+=log(double(m));
 return ret;
 }


//---------------------------------------------------------------------------[ binomial ]-
double log_binomial(unsigned D, unsigned n)
 {
 CHECK_COMPARE( D, >=, n );
 return( log_factorial(D) - log_factorial(n) - log_factorial(D-n) );
 }


//-----------------------------------------------------------------------[ NoFockStates ]-
size_t NoFockStates(size_t N, size_t NoModes, bool PauliF)
 {
 if( !NoModes )  return 0;
 if( !PauliF )   return binomial( NoModes+N-1, NoModes-1 ); 
 CHECK_COMPARE( N, <=, NoModes );
 return binomial( NoModes, N );
 }



//##################################################################################################

//-----------------------------------------------------------------[ operatorBasis_full ]-
Monomials operatorBasis_full(Subsystems const& Omega, int algebraF, sdp_gs::Indices& NoElements)
 {
 DEBUG_RANGE_CHECK( algebraF, 0, 2 );
 bool PauliF = algebraF == 1 || algebraF ==2;
 Monomials A;
 NoElements = sdp_gs::Indices( Omega.NoSubsystems()+1 );
 
 A << Monomial();
 NoElements[0] = 1;
 
 for(size_t k = 0; k < Omega.NoSubsystems(); k++)
  {
  
  for(size_t m = 0; m <= k+1; m++)
   {
   if(k==0)
    { TRACE(m)(k+1-m);  }
   Monomials Ak(   NoFockStates( m,     Omega[k].size(), PauliF )
                 * NoFockStates( k+1-m, Omega[k].size(), PauliF ) );
   size_t cnt = 0;
   for(FockState sm( m, Omega[k].size(), PauliF ); sm.is_valid(); sm.inc())
    for(FockState sn( k+1-m, Omega[k].size(), PauliF ); sn.is_valid(); sn.inc(), cnt++)
     Ak[cnt] = annihilatorProduct( Omega[k][sm()] ) * herm( annihilatorProduct( Omega[k][sn()] ) );
   CHECK_EQUAL( cnt, Ak.size() );
   A = direct_sum( A, Ak );
   }
  
  NoElements[k+1] = A.size();
  }
 
 return A;
 }


//##################################################################################################
