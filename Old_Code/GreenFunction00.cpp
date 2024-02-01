#include "GreenFunction00.h"

#include "latticeAlgebra00.h"
#include "operatorTransform00.h"

#include "linearalgebra/vector_utility.h"	// random_vector



//__________________________________________________________________________________________________
//##################################################################################### GreenBasis #

//-----------------------------------------------------------------------[ stream stuff ]-
template<>
const unsigned sdp_gs::ClassID<GreenBasis>::ID = 7553801;

std::ostream& operator<<(std::ostream& out, GreenBasis const& v)
 {
 out << "(algebraF="     << v.algebraF()     << ",\n"
     << " elementsReal=" << v.elementsReal() << ",\n"
     << " elementsImag=" << v.elementsImag() << ",\n"
     << " repMap.size()=" << v.representativeMap().size()
     << ")";
 return out;
 }

PStream::opstream& GreenBasis::write(PStream::opstream &out) const
 {
 out << (unsigned) sdp_gs::ClassID<GreenBasis>::ID;	// ID
 out << (int)0; 					// version
 out << AlgebraF << Basis << symm << ElementsReal << ElementsImag << repMap;
 return out;
 }

PStream::ipstream& GreenBasis::read(PStream::ipstream &in)
 {
 unsigned ID;
 in >> ID;
 CHECK_EQUAL( ID, sdp_gs::ClassID<GreenBasis>::ID );
 int ver;
 in >> ver;
 DEBUG_CHECK_EQUAL(ver,0);
 in >> AlgebraF >> Basis >> symm >> ElementsReal >> ElementsImag >> repMap;
 return in;
 }


//------------------------------------------------------------------------[ constructor ]-
GreenBasis::GreenBasis(int i_algebraF, OneParticleBasis const& i_basis, Symmetries const& i_symm,
                       MonomialQnList const& i_elementsReal, MonomialQnList const& i_elementsImag,
                       RepresentativeMap const& i_repMap)
  : AlgebraF(i_algebraF), Basis(i_basis), symm(i_symm),
    ElementsReal(i_elementsReal), ElementsImag(i_elementsImag), repMap(i_repMap)
 {
 DEBUG_RANGE_CHECK( AlgebraF, 0, 2 );
 }


//-------------------------------------------------------------------[ collect_elements ]-
// If necessary, add new elements to the Green's function basis, that would be necessary
//  to evaluate the given Polynomial.
void GreenBasis::collect_elements(Polynomial const& i_p, bool use_repMap)
 { collect_elements( polynomialQn(i_p,Basis), use_repMap ); }

void GreenBasis::collect_elements(PolynomialQn const& i_p, bool use_repMap)
 {
 PolynomialQn p;
 if( use_repMap )  p = representative( erased_zeros(i_p,AlgebraF,symm), symm, AlgebraF, repMap );
 else              p = representative( erased_zeros(i_p,AlgebraF,symm), symm, AlgebraF );
 
 for(PolynomialQn::coeff_cit i = p.begin(); i != p.end(); ++i)
  if( std::abs(i->second) > zeroThresh )
   {
   // Check whether Re(G) and/or Im(G) are needed for the evaluation of the polynomial.
   bool need_Re = ( std::abs( p[ herm(i->first) ] + i->second ) > zeroThresh );
   bool need_Im = ( !symm.is_real() &&
                    std::abs( p[ herm(i->first) ] - i->second ) > zeroThresh );
   
   // herm(op) < op ?
   Representative rep, repH;
   if( use_repMap )
    {
    rep  = representative( i->first,       symm, AlgebraF, repMap );
    repH = representative( herm(i->first), symm, AlgebraF, repMap );
    }
   else
    {
    rep  = representative( i->first,       symm, AlgebraF );
    repH = representative( herm(i->first), symm, AlgebraF );
    }
   if( repH() < rep() )  rep = repH;
   
   // Get symmetry orbit to determine symmetry constraints on Re(G) and Im(G).
   //SymmetryOrbit orbit( op, symm, AlgebraF );
   if( rep.is_zero_Re() )  need_Re = 0;
   if( rep.is_zero_Im() )  need_Im = 0;
   
   // Add op to basis, if needed.
   if( need_Re && !ElementsReal.contains( rep() ) )
    ElementsReal.push_back( rep() );
   if( need_Im && !ElementsImag.contains( rep() ) )
    ElementsImag.push_back( rep() );
   }
 }


//------------------------------------------------------------------------[ expand_dual ]-
// G[H] = \sum_s H_s G_s = \sum_s H_s [ Re(G_s) + i*Im(G_s) ]
// p = a*s + b*s^+  =>  G[p] = a*g + b*g^* = Re(g)*(a+b) + Im(g)*i*(a-b)
GreenFunctionDual GreenBasis::expand_dual(Polynomial const& p) const
 { return expand_dual( polynomialQn(p,Basis) ); }

GreenFunctionDual GreenBasis::expand_dual(PolynomialQn const& i_p) const
 {
 GreenFunctionDual ret( sdp_gs::PointCF(sizeReal(),0), sdp_gs::PointCF(sizeImag(),0) );
 PolynomialQn p = representative( erased_zeros(i_p,AlgebraF,symm), symm, AlgebraF, repMap );
 
 for(PolynomialQn::coeff_cit i = p.begin(); i != p.end(); ++i)
  if( std::abs(i->second) > zeroThresh )
   {
   // Check whether Re(G) and/or Im(G) are needed for the evaluation of the polynomial.
   bool need_Re = ( std::abs( p[ herm(i->first) ] + i->second ) > zeroThresh );
   bool need_Im = ( !symm.is_real() &&
                    std::abs( p[ herm(i->first) ] - i->second ) > zeroThresh );
   
   // herm(op) < op ?
   bool use_adjoint = 0;
   Representative rep  = representative( i->first,       symm, AlgebraF, repMap );
   Representative repH = representative( herm(i->first), symm, AlgebraF, repMap );
   if( repH() < rep() )  { rep = repH;  use_adjoint = 1; }
   
   // Get symmetry orbit to determine symmetry constraints on Re(G) and Im(G).
   if( rep.is_zero_Re() )  need_Re = 0;
   if( rep.is_zero_Im() )  need_Im = 0;

   if( need_Re )
    {
    size_t closest;
    CHECK_COMPARE( distance_min( rep(), ElementsReal.labels(), closest ), <, 100*zeroThresh );
    ret.first[closest] += i->second;
    }
   if( need_Im )
    {
    size_t closest;
    CHECK_COMPARE( distance_min( rep(), ElementsImag.labels(), closest ), <, 100*zeroThresh );
    if( !use_adjoint )  ret.second[closest] += i->second;
    else                ret.second[closest] -= i->second;
    }
   }
 return ret;
 }


//-------------------------------------------------------------------------------[ eval ]-
sdp_gs::CFloat GreenBasis::eval(Polynomial const& p, GreenFunction const& G) const
 { return eval( polynomialQn(p,Basis), G ); }

sdp_gs::CFloat GreenBasis::eval(PolynomialQn const& i_p, GreenFunction const& G) const
 {
 sdp_gs::CFloat ret = 0;
 PolynomialQn p = representative( erased_zeros(i_p,AlgebraF,symm), symm, AlgebraF, repMap );
 
 for(PolynomialQn::coeff_cit i = p.begin(); i != p.end(); ++i)
  if( std::abs(i->second) > zeroThresh )
   {
   // Check whether Re(G) and/or Im(G) are needed for the evaluation of the polynomial.
   bool need_Re = ( std::abs( p[ herm(i->first) ] + i->second ) > zeroThresh );
   bool need_Im = ( !symm.is_real() &&
                    std::abs( p[ herm(i->first) ] - i->second ) > zeroThresh );
   
   // herm(op) < op ?
   bool use_adjoint = 0;
   Representative rep  = representative( i->first,       symm, AlgebraF, repMap );
   Representative repH = representative( herm(i->first), symm, AlgebraF, repMap );
   if( repH() < rep() )  { rep = repH;  use_adjoint = 1; }
   
   // Get symmetry orbit to determine symmetry constraints on Re(G) and Im(G).
   if( rep.is_zero_Re() )  need_Re = 0;
   if( rep.is_zero_Im() )  need_Im = 0;

   if( need_Re )
    {
    size_t closest;
    CHECK_COMPARE( distance_min( rep(), ElementsReal.labels(), closest ), <, 100*zeroThresh );
    ret += G.first[closest] * i->second;
    }
   if( need_Im )
    {
    size_t closest;
    CHECK_COMPARE( distance_min( rep(), ElementsImag.labels(), closest ), <, 100*zeroThresh );
    if( !use_adjoint )  ret += sdp_gs::my_i * G.second[closest] * i->second;
    else                ret -= sdp_gs::my_i * G.second[closest] * i->second;
    }
   }
 return ret;
 }


//-------------------------------------------------------------------------------[ eval ]-
sdp_gs::Float GreenBasis::eval_Re(Polynomial const& p, GreenFunction const& G) const
 { return eval_Re( polynomialQn(p,Basis), G ); }

sdp_gs::Float GreenBasis::eval_Re(PolynomialQn const& i_p, GreenFunction const& G) const
 {
 using namespace LinearAlgebra;
 sdp_gs::Float ret = 0;
 PolynomialQn p = representative( erased_zeros(i_p,AlgebraF,symm), symm, AlgebraF, repMap );
 
 for(PolynomialQn::coeff_cit i = p.begin(); i != p.end(); ++i)
  if( std::abs(i->second) > zeroThresh )
   {
   // Check whether Re(G) and/or Im(G) are needed for the evaluation of the polynomial.
   bool need_Re = ( std::abs( p[ herm(i->first) ] + i->second ) > zeroThresh );
   bool need_Im = ( !symm.is_real() &&
                    std::abs( p[ herm(i->first) ] - i->second ) > zeroThresh );
   
   // herm(op) < op ?
   bool use_adjoint = 0;
   Representative rep  = representative( i->first,       symm, AlgebraF, repMap );
   Representative repH = representative( herm(i->first), symm, AlgebraF, repMap );
   if( repH() < rep() )  { rep = repH;  use_adjoint = 1; }
   
   // Get symmetry orbit to determine symmetry constraints on Re(G) and Im(G).
   if( rep.is_zero_Re() )  need_Re = 0;
   if( rep.is_zero_Im() )  need_Im = 0;

   if( need_Re )
    {
    size_t closest;
    CHECK_COMPARE( distance_min( rep(), ElementsReal.labels(), closest ), <, 100*zeroThresh );
    ret += G.first[closest] * sdp_gs::real(i->second);
    }
   if( need_Im )
    {
    size_t closest;
    CHECK_COMPARE( distance_min( rep(), ElementsImag.labels(), closest ), <, 100*zeroThresh );
    if( !use_adjoint )  ret += -G.second[closest] * sdp_gs::imag(i->second);
    else                ret +=  G.second[closest] * sdp_gs::imag(i->second);
    }
   }
 return ret;
 }



//__________________________________________________________________________________________________
//########################################################################### random_GreenFunction #
GreenFunction random_GreenFunction_Re(GreenBasis const& G_basis)
 {
 GreenFunction G;
 G.first  = LinearAlgebra::random_vector<sdp_gs::Float>( G_basis.sizeReal() )
           -sdp_gs::PointF( G_basis.sizeReal(), 0.5 );
 G.second = sdp_gs::PointF( G_basis.sizeImag(), 0 );
 return G;
 }

GreenFunction random_GreenFunction(GreenBasis const& G_basis)
 {
 GreenFunction G;
 G.first  = LinearAlgebra::random_vector<sdp_gs::Float>( G_basis.sizeReal() )
           -sdp_gs::PointF( G_basis.sizeReal(), 0.5 );
 G.second = LinearAlgebra::random_vector<sdp_gs::Float>( G_basis.sizeImag() )
           -sdp_gs::PointF( G_basis.sizeImag(), 0.5 );
 return G;
 }



//##################################################################################################
