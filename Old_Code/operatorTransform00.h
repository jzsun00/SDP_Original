//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 15.07.2011 Thomas Barthel
//% SDP groundstate kit
//% Transformation of quantum numbers and operators under symmetry operations.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// key words: TODO, QUESTION, TEST, CHECK, REMARK, CAUTION
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#ifndef SDP_GS_KIT_OPERATORTRANSFORM00_H
#define SDP_GS_KIT_OPERATORTRANSFORM00_H

#include "algebra/operator00.h"
#include "lattice/symmetries00.h"



//__________________________________________________________________________________________________
//############################################################### transform operator to operatorQn #

//-------------------------------------------------------------------------[ ladderOpQn ]-
inline LadderOpQn ladderOpQn(LadderOp const& op, OneParticleBasis const& basis)
 { return LadderOpQn( basis[op.index()], op.is_creator() ); }


//-------------------------------------------------------------------------[ monomialQn ]-
MonomialQn monomialQn(Monomial const& op, OneParticleBasis const& basis);


//------------------------------------------------------------------------[ monomialsQn ]-
template <typename T>
MonomialsQn monomialsQn(T const& op, OneParticleBasis const& basis,
                        typename boost::enable_if< LinearAlgebra::is_vector<T> >::type* dummy=0)
 {
 MonomialsQn ret(op.size());
 for(size_t n = 0; n < op.size(); n++)  ret[n] = monomialQn( Monomial(op[n]), basis );
 return ret;
 }

inline MonomialsQn const& monomialsQn(MonomialsQn const& op, OneParticleBasis const& basis)
 { return op; }


//-----------------------------------------------------------------------[ polynomialQn ]-
PolynomialQn polynomialQn(Polynomial const& op, OneParticleBasis const& basis);

inline PolynomialQn const& polynomialQn(PolynomialQn const& op, OneParticleBasis const& basis)
 { return op; }


//----------------------------------------------------------------------[ polynomialsQn ]-
template <typename T>
PolynomialsQn polynomialsQn(T const& op, OneParticleBasis const& basis,
                            typename boost::enable_if< LinearAlgebra::is_vector<T> >::type* dummy=0)
 {
 PolynomialsQn ret(op.size());
 for(size_t n = 0; n < op.size(); n++)  ret[n] = polynomialQn( Polynomial(op[n]), basis );
 return ret;
 }

inline PolynomialsQn const& polynomialsQn(PolynomialsQn const& op, OneParticleBasis const& basis)
 { return op; }



//__________________________________________________________________________________________________
//############################################################### transform operatorQn to operator #

//-----------------------------------------------------------[ collect_oneParticleBasis ]-
inline void collect_oneParticleBasis(OneParticleBasis& basis, LadderOpQn const& op)
 {  basis.insert( op.index() ); }
inline void collect_oneParticleBasis(OneParticleBasis& basis, MonomialQn const& op)
 { for(size_t n = 0; n < op.size(); n++)  collect_oneParticleBasis( basis, op[n] ); }

inline void collect_oneParticleBasis(OneParticleBasis& basis, PolynomialQn const& op)
 { for(PolynomialQn::coeff_cit i = op.begin(); i != op.end(); ++i)
    collect_oneParticleBasis( basis, i->first );  }

inline void collect_oneParticleBasis(OneParticleBasis& basis, PolynomialsQn const& op)
 { for(size_t n = 0; n < op.size(); n++)  collect_oneParticleBasis( basis, op[n] ); }


//---------------------------------------------------------------------------[ ladderOp ]-
inline LadderOp ladderOp(LadderOpQn const& op, OneParticleBasis const& basis)
 { return LadderOp( basis(op.index()), op.is_creator() ); }


//---------------------------------------------------------------------------[ monomial ]-
Monomial monomial(MonomialQn const& op, OneParticleBasis const& basis);


//-------------------------------------------------------------------------[ polynomial ]-
Polynomial polynomial(PolynomialQn const& op, OneParticleBasis const& basis);

inline Polynomial const& polynomial(Polynomial const& op, OneParticleBasis const& basis)
 { return op; }


//------------------------------------------------------------------------[ polynomials ]-
template <typename T>
Polynomials polynomials(T const& op, OneParticleBasis const& basis,
                        typename boost::enable_if< LinearAlgebra::is_vector<T> >::type* dummy=0)
 {
 Polynomials ret(op.size());
 for(size_t n = 0; n < op.size(); n++)  ret[n] = polynomial( PolynomialQn(op[n]), basis );
 return ret;
 }

inline Polynomials const& polynomials(Polynomials const& op, OneParticleBasis const& basis)
 { return op; }



//__________________________________________________________________________________________________
//########################################################## finite symmetry group transformations #

//--------------------------------------------------------------------------[ transform ]-
inline void transform(QuantumNumber& qn, sdp_gs::Indices const& trafo_idcs,
                      Symmetries::operator_type const& trafo)
 { DEBUG_CHECK_EQUAL( trafo.size1(), trafo_idcs.size() );
   qn[trafo_idcs] = trafo*qn[trafo_idcs]; }

inline QuantumNumber transformed(QuantumNumber const& qn, sdp_gs::Indices const& trafo_idcs,
                                 Symmetries::operator_type const& trafo)
 { QuantumNumber ret(qn);  transform(ret,trafo_idcs,trafo);  return ret; }


inline void transform(LadderOpQn& op, sdp_gs::Indices const& trafo_idcs,
                            Symmetries::operator_type const& trafo)
 { transform(op.index(),trafo_idcs,trafo); }

inline LadderOpQn transformed(LadderOpQn const& op, sdp_gs::Indices const& trafo_idcs,
                              Symmetries::operator_type const& trafo)
 { LadderOpQn ret(op);  transform(ret,trafo_idcs,trafo);  return ret; }


inline void transform(MonomialQn& op, sdp_gs::Indices const& trafo_idcs,
                      Symmetries::operator_type const& trafo)
 { for(size_t n = 0; n < op.size(); n++)  op[n] = transformed(op[n],trafo_idcs,trafo); }

inline MonomialQn transformed(MonomialQn const& op, sdp_gs::Indices const& trafo_idcs,
                              Symmetries::operator_type const& trafo)
 { MonomialQn ret(op);  transform(ret,trafo_idcs,trafo);  return ret; }





//__________________________________________________________________________________________________
//########################################################################### lattice translations #

//--------------------------------------------------------------------------[ translate ]-
inline void translate(QuantumNumber& qn, sdp_gs::Indices const& coord_idcs,
                      Lattice const& latt, sdp_gs::PointI const& shift)
 { qn[coord_idcs] += latt.basis()*shift; }

inline QuantumNumber translated(QuantumNumber const& qn, sdp_gs::Indices const& coord_idcs,
                                Lattice const& latt, sdp_gs::PointI const& shift)
 { QuantumNumber ret(qn);  translate(ret,coord_idcs,latt,shift);  return ret; }


inline void translate(LadderOpQn& op, sdp_gs::Indices const& coord_idcs,
                      Lattice const& latt, sdp_gs::PointI const& shift)
 { translate(op.index(),coord_idcs,latt,shift); }

inline LadderOpQn translated(LadderOpQn const& op, sdp_gs::Indices const& coord_idcs,
                             Lattice const& latt, sdp_gs::PointI const& shift)
 { LadderOpQn ret(op);  translate(ret,coord_idcs,latt,shift);  return ret; }


inline void translate(MonomialQn& op, sdp_gs::Indices const& coord_idcs,
                      Lattice const& latt, sdp_gs::PointI const& shift)
 { for(size_t n = 0; n < op.size(); n++)  translate(op[n],coord_idcs,latt,shift); }

inline MonomialQn translated(MonomialQn const& op, sdp_gs::Indices const& coord_idcs,
                             Lattice const& latt, sdp_gs::PointI const& shift)
 { MonomialQn ret(op);  translate(ret,coord_idcs,latt,shift);  return ret; }


PolynomialQn translated(PolynomialQn const& op, sdp_gs::Indices const& coord_idcs,
                        Lattice const& latt, sdp_gs::PointI const& shift);


//---------------------------------------------------------------------[ standard_shift ]-
sdp_gs::PointI standard_shift(QuantumNumber const& qn, sdp_gs::Indices const& coord_idcs,
                              Lattice const& latt);

inline sdp_gs::PointI standard_shift(LadderOpQn const& op, sdp_gs::Indices const& coord_idcs,
                                     Lattice const& latt)
 { return standard_shift( op.index(), coord_idcs, latt ); }

inline sdp_gs::PointI standard_shift(MonomialQn const& op, sdp_gs::Indices const& coord_idcs,
                                     Lattice const& latt)
 { if( !op.size() )  return( sdp_gs::PointI( latt.basis().size2(), 0 ) );
   return standard_shift( op[0], coord_idcs, latt ); }


//----------------------------------------------------------------------[ translate_std ]-
inline void translate_std(MonomialQn& op, sdp_gs::Indices const& coord_idcs, Lattice const& latt)
 { translate( op, coord_idcs, latt, standard_shift( op, coord_idcs, latt ) ); }

inline MonomialQn translated_std(MonomialQn const& op, sdp_gs::Indices const& coord_idcs, Lattice const& latt)
 { MonomialQn ret(op);  translate_std(ret,coord_idcs,latt);  return ret; }



//__________________________________________________________________________________________________
//################################################################################## SymmetryOrbit #
class SymmetryOrbit
 {
 public:
 SymmetryOrbit(MonomialQn const& i_op, Symmetries const& symm, int algebraF);
 
 size_t size() const  { return op.size(); }
 
 MonomialsQn    const& operator()()           const  { return op; }
 MonomialQn     const& operator[](size_t n)   const  { return op[n]; }
 sdp_gs::PointI const& coeffs()               const  { return Coeff; }
 int                   coeff     (size_t n)   const  { return Coeff[n]; }
 int                   representative_coeff() const  { return Coeff[nRep]; }
 
 MonomialQn     const& representative()     const  { return op[nRep]; }
 int                   zeroF()              const  { return ZeroF; }
 bool                  is_zero()            const  { return( (ZeroF&3) == 3 ); }
 bool                  is_zero_Re()         const  { return( (ZeroF&1) == 1 ); }
 bool                  is_zero_Im()         const  { return( (ZeroF&2) == 2 ); }
 
 MonomialsQn           alternatives() const;  // all operators from the orbit that are different from the rep.
 
 protected:
 MonomialsQn    op;
 sdp_gs::PointI Coeff;
 size_t         nRep;
 int            ZeroF;	// bit1: Re(G[op])=0,  bit2: Im(G[op])=0
 };



//__________________________________________________________________________________________________
//################################################################################# Representative #
// This is a, hopefully, more efficient version of SymmetryOrbit.
// Different from SymmetryOrbit, Representative does not give access to all elements of the orbit,
//  but only the minimal element, the representative.
class Representative
 {
 public:
 Representative()  {}
 Representative(MonomialQn const& i_op, int i_coeff, int i_zeroF)
  : op(i_op), Coeff(i_coeff), ZeroF(i_zeroF)  {}
 Representative& operator=(Representative const& v)
  { op = v.op;  Coeff = v.Coeff;  ZeroF = v.ZeroF;  return *this; }
 
 MonomialQn const& operator()() const  { return op; }
 int               coeff()      const  { return Coeff; }
 int               zeroF()      const  { return ZeroF; }
 bool              is_zero()    const  { return( (ZeroF&3) == 3 ); }
 bool              is_zero_Re() const  { return( (ZeroF&1) == 1 ); }
 bool              is_zero_Im() const  { return( (ZeroF&2) == 2 ); }
 
 void set_coeff(int i_coeff)  { Coeff = i_coeff; }
 
 PStream::opstream& write(PStream::opstream &out) const  { out << op << Coeff << ZeroF;  return out; }
 PStream::ipstream& read (PStream::ipstream &in)         { in  >> op >> Coeff >> ZeroF;  return in; }
 
 protected:
 MonomialQn op;
 int        Coeff;
 int        ZeroF;	// bit1: Re(G[op])=0,  bit2: Im(G[op])=0
 };


//-----------------------------------------------------------------------[ stream stuff ]-
inline PStream::opstream& operator<<(PStream::opstream& out, Representative const& v)
 { return v.write(out); }
inline PStream::ipstream& operator>>(PStream::ipstream& in , Representative      & v)
 { return v.read (in);  }



//__________________________________________________________________________________________________
//############################################################################## RepresentativeMap #
typedef std::map<MonomialQn,Representative> RepresentativeMap;




//__________________________________________________________________________________________________
//################################################################### QuantumNumber representative #

//---------------------------------------------------------------------[ representative ]-
QuantumNumber representative(QuantumNumber const& qn, Symmetries const& symm);



//__________________________________________________________________________________________________
//######################################################################## monomial representative #

//---------------------------------------------------------------------[ representative ]-
Representative representative(MonomialQn const& op, Symmetries const& symm, int algebraF,
                              RepresentativeMap const& repMap);

Representative representative(MonomialQn const& op, Symmetries const& symm, int algebraF,
                              RepresentativeMap& repMap);

inline Representative representative(MonomialQn const& op, Symmetries const& symm, int algebraF)
 { return representative( op, symm, algebraF, RepresentativeMap() ); }



//__________________________________________________________________________________________________
//###################################################################### polynomial representative #

//---------------------------------------------------------------------[ representative ]-
// Replace all Monomials by a representative in a standardized fashion.
PolynomialQn representative(PolynomialQn const& op, Symmetries const& symm, int algebraF,
                            RepresentativeMap const& repMap);

inline PolynomialQn representative(PolynomialQn const& op, Symmetries const& symm, int algebraF)
 { return representative( op, symm, algebraF, RepresentativeMap() ); }



//##################################################################################################

#endif // SDP_GS_KIT_OPERATORTRANSFORM00_H
