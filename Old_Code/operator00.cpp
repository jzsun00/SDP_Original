#include "operator00.h"



//__________________________________________________________________________________________________
//####################################################################################### Monomial #

//----------------------------------------------------------------------[ quantumNumber ]-
QuantumNumber quantumNumber(MonomialQn const& x, sdp_gs::Indices const& qn_idcs)
 {
 LadderOpQn::index_type qn( qn_idcs.size(), 0 );
 for(size_t n = 0; n < x.size(); n++)
  {
  if( x[n].is_creator() )  qn += x[n].index()[qn_idcs];
  else                     qn -= x[n].index()[qn_idcs];
  }
 return qn;
 }



//__________________________________________________________________________________________________
//##################################################################################### Polynomial #

//----------------------------------------------------------------------[ quantumNumber ]-
QuantumNumber quantumNumber(PolynomialQn const& x, sdp_gs::Indices const& qn_idcs)
 {
 if( !x.size() )  return QuantumNumber( qn_idcs.size(), 0 );
 PolynomialQn::coeff_cit i = x.begin();
 QuantumNumber qn = quantumNumber( i->first, qn_idcs );
 for(++i; i != x.end(); ++i)
  if( qn != quantumNumber( i->first, qn_idcs ) )  return QuantumNumber();
 return qn;
 }



//__________________________________________________________________________________________________
//######################################################################################### random #
Monomial random_Monomial(size_t NoModes, size_t degree)
 {
 Monomial ret;
 for(size_t d = 0; d < degree; d++)  ret *= random_LadderOp( NoModes );
 return ret;
 }


Polynomial random_Polynomial(size_t NoModes, sdp_gs::Indices const& NoTerms)
 {
 Polynomial ret;
 ret += ( double(rand())/RAND_MAX-0.5 + sdp_gs::my_i*(double(rand())/RAND_MAX-0.5) ) * Monomial();
 for(size_t d = 0; d < NoTerms.size(); d++)
  ret += ( double(rand())/RAND_MAX-0.5 + sdp_gs::my_i*(double(rand())/RAND_MAX-0.5) )*random_Monomial( NoModes, d );
 return ret;
 }


Monomial random_Monomial(sdp_gs::Indices const& subsystem, size_t degree)
 {
 Monomial ret;
 for(size_t d = 0; d < degree; d++)  ret *= random_LadderOp( subsystem );
 return ret;
 }


Polynomial random_Polynomial(LinearAlgebra::Vector<sdp_gs::Indices> const& subsystems,
                             sdp_gs::Indices const& NoTerms)
 {
 CHECK_EQUAL( subsystems.size(), NoTerms.size() );
 Polynomial ret;
 ret += ( double(rand())/RAND_MAX-0.5 + sdp_gs::my_i*(double(rand())/RAND_MAX-0.5) ) * Monomial();
 for(size_t d = 0; d < NoTerms.size(); d++)
  for(size_t nt = 0; nt < NoTerms[d]; nt++)
   ret += ( double(rand())/RAND_MAX-0.5 + sdp_gs::my_i*(double(rand())/RAND_MAX-0.5) )
          * random_Monomial( subsystems[d], d+1 );
 return ret;
 }
