//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 29.06.2011 Thomas Barthel
//% SDP groundstate kit
//% Evaluation of Monomials and Polynomials with respect to Green's functions.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// key words: TODO, QUESTION, TEST, CHECK, REMARK, CAUTION
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#ifndef SDP_GS_KIT_GREENFUNCTION00_H
#define SDP_GS_KIT_GREENFUNCTION00_H

// #include "lattice/symmetries00.h"
// #include "algebra/operator00.h"
#include "operatorTransform00.h"



//__________________________________________________________________________________________________
//################################################################################## GreenFunction #
// Note, in principle, we allow here for freely choosing G[Id]. In most cases, one might want G[Id]=1.
typedef std::pair<sdp_gs::PointF, sdp_gs::PointF>  GreenFunction;
typedef std::pair<sdp_gs::PointCF,sdp_gs::PointCF> GreenFunctionDual;
typedef std::pair<sdp_gs::PointF, sdp_gs::PointF>  GreenFunctionDualHerm;


//-----------------------------------------------------------------------[ stream stuff ]-
template <typename T>
std::ostream& operator<<(std::ostream& out, std::pair<T,T> const& v)
 { out << "(" << v.first << ", " << v.second << ")";  return out; }


//-------------------------------------------------------------------------[ inner_prod ]-
inline sdp_gs::CFloat inner_prod(GreenFunctionDual const& p, GreenFunction const& G)
 { //TRACE(LinearAlgebra::parallel_prod( p.first,  G.first ))(sdp_gs::my_i*LinearAlgebra::parallel_prod( p.second, G.second ))(p.second[LinearAlgebra::range(0,8)])(G.second[LinearAlgebra::range(0,8)]);
  return                LinearAlgebra::parallel_prod( p.first,  G.first ) 
         + sdp_gs::my_i*LinearAlgebra::parallel_prod( p.second, G.second ); }

//----------------------------------------------------------------------[ inner_prod_Re ]-
inline sdp_gs::Float inner_prod_Re(GreenFunctionDual const& p, GreenFunction const& G)
 { return   LinearAlgebra::parallel_prod( LinearAlgebra::real(p.first),  G.first )
          - LinearAlgebra::parallel_prod( LinearAlgebra::imag(p.second), G.second ); }

//--------------------------------------------------------------[ greenFunctionDualHerm ]-
inline GreenFunctionDualHerm greenFunctionDualHerm(GreenFunctionDual const& p)
 { DEBUG_CHECK_COMPARE( normInf(imag(p.first)),  <, zeroThresh );
   DEBUG_CHECK_COMPARE( normInf(real(p.second)), <, zeroThresh );
   return GreenFunctionDualHerm( real(p.first), -imag(p.second) ); }


//__________________________________________________________________________________________________
//##################################################################################### GreenBasis #
// Basis for the (required) Green's functions.
// One "symmetry" that is always given, is Hermitian conjugation, i.e., we need only one basis element,
//  either \sigma or \sigma^+.
class GreenBasis
 {
 public:
 typedef LabelList<MonomialQn> MonomialQnList;
 typedef MonomialQnList::map_type map_type;
 
 GreenBasis() : AlgebraF(-1)  {}
 GreenBasis(int i_algebraF, OneParticleBasis const& i_basis, Symmetries const& i_symm)
  : AlgebraF(i_algebraF), Basis(i_basis), symm(i_symm)  {}
 GreenBasis(int i_algebraF, OneParticleBasis const& i_basis, Symmetries const& i_symm,
            MonomialQnList const& i_elementsReal, MonomialQnList const& i_elementsImag,
            RepresentativeMap const& i_repMap);
 GreenBasis& operator=(GreenBasis const& v)
  { AlgebraF = v.AlgebraF;  Basis = v.Basis;  symm = v.symm;
    ElementsReal = v.ElementsReal;  ElementsImag = v.ElementsImag;  repMap = v.repMap;  return *this; }
 bool operator==(GreenBasis const& v) const  // "!(==)" is not (yet) "!="
  { return( AlgebraF == v.AlgebraF && symm == v.symm
            && ElementsReal == v.ElementsReal && ElementsImag == v.ElementsImag ); }
 
 size_t sizeReal() const  { return ElementsReal.size(); }
 size_t sizeImag() const  { return ElementsImag.size(); }
 size_t size()     const  { return( sizeReal() + sizeImag() ); }
 
 int                     algebraF()       const  { return AlgebraF; }
 Symmetries       const& symmetries()     const  { return symm; }
 OneParticleBasis const& basis()          const  { return Basis; }
 MonomialQnList     const& elementsReal() const  { return ElementsReal; }
 MonomialQnList     const& elementsImag() const  { return ElementsImag; }
 RepresentativeMap  const& representativeMap() const  { return repMap; }
 
 // The following two functions can be used to probe for cases where basis elements occur
 //  multiple times due to small numerical deviations in the associated quantum numbers.
 double distance_min_Re() const  { return distance_min( ElementsReal.labels() ); }
 double distance_min_Im() const  { return distance_min( ElementsImag.labels() ); }
 
 void              collect_elements(PolynomialQn const& p, bool use_repMap=1);
 GreenFunctionDual expand_dual(PolynomialQn const& p) const;
 sdp_gs::CFloat    eval       (PolynomialQn const& p, GreenFunction const& G) const;
 sdp_gs::Float     eval_Re    (PolynomialQn const& p, GreenFunction const& G) const;
 
 // In the following variants, the Polynomials must be given in the single-particle basis Basis.
 void              collect_elements(Polynomial const& p, bool use_repMap=1);
 GreenFunctionDual expand_dual(Polynomial const& p) const;
 sdp_gs::CFloat    eval       (Polynomial const& p, GreenFunction const& G) const;
 sdp_gs::Float     eval_Re    (Polynomial const& p, GreenFunction const& G) const;
//  sdp_gs::Float  eval_Re(Polynomial const& p, GreenFunction const& G) const
//   { return real(eval(p,G)); }
 
 PStream::opstream& write(PStream::opstream &out) const;
 PStream::ipstream& read (PStream::ipstream &in);
 
 protected:
 int              AlgebraF;
 OneParticleBasis Basis;
 Symmetries       symm;
 
 MonomialQnList ElementsReal;	// representatives of the required MonomialQn's
 MonomialQnList ElementsImag;
 
 RepresentativeMap repMap;
 };


//-----------------------------------------------------------------------[ stream stuff ]-
inline PStream::opstream& operator<<(PStream::opstream& out, GreenBasis const& v)
 { return v.write(out); }
inline PStream::ipstream& operator>>(PStream::ipstream& in , GreenBasis      & v)
 { return v.read (in);  }
std::ostream& operator<<(std::ostream& out, GreenBasis const& v);



//__________________________________________________________________________________________________
//########################################################################### random_GreenFunction #
GreenFunction random_GreenFunction_Re(GreenBasis const& G_basis);
GreenFunction random_GreenFunction   (GreenBasis const& G_basis);



//##################################################################################################

#endif // SDP_GS_KIT_GREENFUNCTION00_H
