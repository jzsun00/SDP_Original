#include "QuantumNumber00.h"


//##################################################################################################

//-----------------------------------------------------------------------[ qn_projected ]-
QuantumNumbers qn_projected(QuantumNumbers const& qns, sdp_gs::Indices const& idcs)
 {
 QuantumNumbers ret( qns );
 for(size_t n = 0; n < ret.size(); n++)
  ret[n] = ret[n][idcs];
 return ret;
 }


//-----------------------------------------------------------------------[ qn_directSum ]-
QuantumNumbers qn_directSum(QuantumNumbers const& qns1, QuantumNumbers const& qns2)
 {
 CHECK_EQUAL( qns1.size(), qns2.size() );
 if( !qns1.size() )  return qns1;
 
 QuantumNumbers ret( qns1.size(), QuantumNumber( qns1[0].size()+qns2[0].size() ) );
 for(size_t n = 0; n < ret.size(); n++)
  ret[n] = LinearAlgebra::direct_sum( qns1[n], qns2[n] );
 return ret;
 }



//##################################################################################################
