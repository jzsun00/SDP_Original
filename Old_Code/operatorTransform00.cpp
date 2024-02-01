#include "operatorTransform00.h"

#include "algebra/algebra00.h"

#include "linalg/TensorL02.h"



//__________________________________________________________________________________________________
//############################################################### transform operator to operatorQn #

//-------------------------------------------------------------------------[ monomialQn ]-
MonomialQn monomialQn(Monomial const& op, OneParticleBasis const& basis)
 {
 MonomialQn::container_type ret( op.size() );
 for(size_t n = 0; n < op.size(); n++)  ret[n] = ladderOpQn( op[n], basis );
 return MonomialQn( ret );
 }


//-----------------------------------------------------------------------[ polynomialQn ]-
PolynomialQn polynomialQn(Polynomial const& op, OneParticleBasis const& basis)
 {
 PolynomialQn ret;
 for(Polynomial::coeff_cit i = op.begin(); i != op.end(); ++i)
  ret[ monomialQn( i->first, basis ) ] = i->second;
 return ret;
 }



//__________________________________________________________________________________________________
//############################################################### transform operatorQn to operator #

//---------------------------------------------------------------------------[ monomial ]-
Monomial monomial(MonomialQn const& op, OneParticleBasis const& basis)
 {
 Monomial::container_type ret( op.size() );
 for(size_t n = 0; n < op.size(); n++)  ret[n] = ladderOp( op[n], basis );
 return Monomial( ret );
 }


//-------------------------------------------------------------------------[ polynomial ]-
Polynomial polynomial(PolynomialQn const& op, OneParticleBasis const& basis)
 {
 Polynomial ret;
 for(PolynomialQn::coeff_cit i = op.begin(); i != op.end(); ++i)
  ret[ monomial( i->first, basis ) ] = i->second;
 return ret;
 }



//__________________________________________________________________________________________________
//########################################################################### lattice translations #

//--------------------------------------------------------------------------[ translate ]-
PolynomialQn translated(PolynomialQn const& op, sdp_gs::Indices const& coord_idcs,
                        Lattice const& latt, sdp_gs::PointI const& shift)
 {
 PolynomialQn ret;
 for(PolynomialQn::coeff_cit i = op.begin(); i != op.end(); ++i)
  ret[ translated( i->first, coord_idcs, latt, shift ) ] = i->second;
 return ret;
 }


//---------------------------------------------------------------------[ standard_shift ]-
sdp_gs::PointI standard_shift(QuantumNumber const& qn, sdp_gs::Indices const& coord_idcs,
                              Lattice const& latt)
 {
 using namespace LinearAlgebra;
 DEBUG_CHECK_EQUAL( latt.basis().size1(), coord_idcs.size() );
 sdp_gs::PointI ret( latt.basis().size2() );
 QuantumNumber pos = qn[coord_idcs];
 for(size_t n = 0; n < ret.size(); n++)
  {
  sdp_gs::Float dist = inner_prod( latt.basis()(all,n), pos ) / norm_frob_sq( latt.basis()(all,n) );
  int shift = -(int)floor( dist );
  ret[n] = shift;
  pos += shift*latt.basis()(all,n);
  }
 return ret;
 }



//__________________________________________________________________________________________________
//################################################################################## SymmetryOrbit #

//------------------------------------------------------------------------[ constructor ]-
SymmetryOrbit::SymmetryOrbit(MonomialQn const& i_op, Symmetries const& symm, int algebraF)
 {
 using namespace LinearAlgebra;
 sdp_gs::Indices const& coordIdcs( symm.lattice().idcs );
 Lattice         const& latt     ( symm.lattice().symm );
 
 nRep  = 0;
 ZeroF = 0;
 
 PRECONDITION( is_normalOrdered( i_op, algebraF ) );
 MonomialQn i_Op  = translated_std( i_op,       coordIdcs, latt );
 MonomialQn i_OpH = translated_std( herm(i_op), coordIdcs, latt );
 
 TensorIndex group_ti( symm.NoGroupElements() );
 op    = MonomialsQn(    symm.NoGroupElementsTot() );
 Coeff = sdp_gs::PointI( symm.NoGroupElementsTot() );
 
 for(TensorIndexState i = group_ti.begin(); i.is_valid(); ++i)
  {
  MonomialQn Op(i_Op);
  for(size_t ng = 0; ng < symm.NoGroups(); ng++)
   transform( Op, symm.group(ng).idcs, symm.group(ng).symm[ i[ng] ] );
  
  PolynomialQn p = normalOrdered( PolynomialQn(Op), algebraF );
  CHECK_EQUAL( p.size(), 1 )( Op )( p );
  CHECK_EQUAL( p.begin()->second, sdp_gs::CFloat(int(sdp_gs::real( p.begin()->second ))) );
  
  size_t cnt = i.count();
  op[cnt]    = p.begin()->first;
  translate_std( op[cnt], coordIdcs, latt );
  Coeff[cnt] = int(sdp_gs::real( p.begin()->second ));
  if( op[cnt] < op[nRep] )  nRep = cnt;
  if( op[cnt] == i_Op && Coeff[cnt] != 1 )  ZeroF |= 1+2;
  if( op[cnt] == i_OpH && Coeff[cnt] != 1  )  ZeroF |= 1;
  if( op[cnt] == i_OpH && Coeff[cnt] != -1 )  ZeroF |= 2;
  }
 }


//-----------------------------------------------------------------------[ alternatives ]-
MonomialsQn SymmetryOrbit::alternatives() const
 {
 using namespace LinearAlgebra;
 MonomialsQn ret( size() );
 size_t cnt = 0;
 for(size_t n = 0; n < size(); n++)
  if(    Distance( representative(), op[n] ) > 100*zeroThresh 
      && distance_min( op[n], ret[range(0,cnt)] ) > 100*zeroThresh )
   {
   ret[cnt] = op[n];
   cnt++;
   }
 return ret[range(0,cnt)];
 }



//__________________________________________________________________________________________________
//################################################################################# Representative #

//---------------------------------------------------------------------[ representative ]-
QuantumNumber representative(QuantumNumber const& i_qn, Symmetries const& symm)
 {
 using namespace LinearAlgebra;
 QuantumNumber ret = i_qn;
 
 TensorIndex group_ti( symm.NoGroupElements() );
 for(TensorIndexState i = group_ti.begin(); i.is_valid(); ++i)
  {
  QuantumNumber qn(i_qn);
  for(size_t ng = 0; ng < symm.NoGroups(); ng++)
   transform( qn, symm.group(ng).idcs, symm.group(ng).symm[ i[ng] ] );
  if( qn < ret )  ret = qn;
  }
 
 return ret;
 }


//-----------------------------------------------------------------[ get_representative ]-
Representative get_representative(MonomialQn const& i_Op, MonomialQn const& i_OpH,
                                  Symmetries const& symm, int algebraF)
 {
 using namespace LinearAlgebra;
 sdp_gs::Indices const& coordIdcs( symm.lattice().idcs );
 Lattice         const& latt     ( symm.lattice().symm );
 
 MonomialQn op = i_Op;
 int Coeff = 1;
 int ZeroF = 0;
 
 TensorIndex group_ti( symm.NoGroupElements() );
 for(TensorIndexState i = group_ti.begin(); i.is_valid(); ++i)
  {
  MonomialQn Op(i_Op);
  for(size_t ng = 0; ng < symm.NoGroups(); ng++)
   transform( Op, symm.group(ng).idcs, symm.group(ng).symm[ i[ng] ] );
  
  PolynomialQn p = normalOrdered( PolynomialQn(Op), algebraF );
  CHECK_EQUAL( p.size(), 1 )( Op )( p );
  CHECK_EQUAL( p.begin()->second, sdp_gs::CFloat(int(sdp_gs::real( p.begin()->second ))) );
  
  Op = translated_std( p.begin()->first, coordIdcs, latt );
  int c = int(sdp_gs::real( p.begin()->second ));
  if( Op < op )
   {
   op = Op;
   Coeff = c;
   }
  if( Op == i_Op  && c !=  1 )  ZeroF |= 1+2;
  if( Op == i_OpH && c !=  1 )  ZeroF |= 1;
  if( Op == i_OpH && c != -1 )  ZeroF |= 2;
  }
 
 return Representative( op, Coeff, ZeroF );
 }


//---------------------------------------------------------------------[ representative ]-
Representative representative(MonomialQn const& i_op, Symmetries const& symm, int algebraF,
                              RepresentativeMap const& repMap)
 {
 sdp_gs::Indices const& coordIdcs( symm.lattice().idcs );
 Lattice         const& latt     ( symm.lattice().symm );
 
 PRECONDITION( is_normalOrdered( i_op, algebraF ) );
 MonomialQn i_Op  = translated_std( i_op,       coordIdcs, latt );
 MonomialQn i_OpH = translated_std( herm(i_op), coordIdcs, latt );
 
 RepresentativeMap::const_iterator i = repMap.find( i_Op );
 if( i != repMap.end() )  return i->second;
 
 return get_representative( i_Op, i_OpH, symm, algebraF );
 }


//---------------------------------------------------------------------[ representative ]-
Representative representative(MonomialQn const& i_op, Symmetries const& symm, int algebraF,
                              RepresentativeMap& repMap)
 {
 using namespace LinearAlgebra;
 sdp_gs::Indices const& coordIdcs( symm.lattice().idcs );
 Lattice         const& latt     ( symm.lattice().symm );
 
 PRECONDITION( is_normalOrdered( i_op, algebraF ) );
 MonomialQn i_Op  = translated_std( i_op,       coordIdcs, latt );
 MonomialQn i_OpH = translated_std( herm(i_op), coordIdcs, latt );
 
 RepresentativeMap::const_iterator i = repMap.find( i_Op );
 if( i != repMap.end() )  return i->second;
 
 Representative rep = get_representative( i_Op, i_OpH, symm, algebraF );
 Representative repC( rep );
 
 TensorIndex group_ti( symm.NoGroupElements() );
 for(TensorIndexState i = group_ti.begin(); i.is_valid(); ++i)
  {
  MonomialQn Op(i_Op);
  for(size_t ng = 0; ng < symm.NoGroups(); ng++)
   transform( Op, symm.group(ng).idcs, symm.group(ng).symm[ i[ng] ] );
  
  PolynomialQn p = normalOrdered( PolynomialQn(Op), algebraF );
  CHECK_EQUAL( p.size(), 1 )( Op )( p );
  CHECK_EQUAL( p.begin()->second, sdp_gs::CFloat(int(sdp_gs::real( p.begin()->second ))) );
  
  Op = translated_std( p.begin()->first, coordIdcs, latt );
  int c = int(sdp_gs::real( p.begin()->second ));
  
  CHECK_EQUAL( abs(c), 1 );
  repC.set_coeff( rep.coeff()/c );
  //repMap[ Op ] = repC;
  repMap.insert( RepresentativeMap::value_type( Op, repC ) );
  }
 
 return rep;
 }



//__________________________________________________________________________________________________
//######################################################################## monomial representative #

//---------------------------------------------------------------------[ representative ]-
PolynomialQn representative(PolynomialQn const& i_p, Symmetries const& symm, int algebraF,
                            RepresentativeMap const& repMap)
 {
 PolynomialQn p = normalOrdered( i_p, algebraF );
 PolynomialQn ret;
 
 for(PolynomialQn::coeff_cit i = p.begin(); i != p.end(); ++i)
  {
  Representative rep  = representative( i->first,       symm, algebraF, repMap );
  Representative repH = representative( herm(i->first), symm, algebraF, repMap );
  // coeff's are integer. Hence, we don't write conj(coeffH).
  if( repH() < rep() )  ret[ herm(repH()) ] += sdp_gs::CFloat(repH.coeff()) * i->second;
  else                  ret[ rep() ]        += sdp_gs::CFloat(rep.coeff())  * i->second;
  }
 return ret;
 }



//##################################################################################################
