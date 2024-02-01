//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 29.06.2011 Thomas Barthel
//% SDP groundstate kit
//% Ladder operators, monomials, and polynomials (in second quantization).
//% [Partly copied from exact/Green-rho_quad.]
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// key words: TODO, QUESTION, TEST, CHECK, REMARK, CAUTION
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#ifndef SDP_GS_KIT_OPERATOR00_H
#define SDP_GS_KIT_OPERATOR00_H

#include "configuration.h"
#include "basis/QuantumNumber00.h"

#include "linearalgebra/pstreamio.h"



//__________________________________________________________________________________________________
//####################################################################################### LadderOp #
template <typename index_type>
class LadderOpT;

typedef LadderOpT<sdp_gs::Index> LadderOp;
typedef LadderOpT<QuantumNumber> LadderOpQn;



//__________________________________________________________________________________________________
//####################################################################################### LadderOp #
template <typename Idx>
class LadderOpT
 {
 public:
 typedef Idx index_type;
 
 LadderOpT() : idx( idx_invalid<index_type>() ), creatorF(0)  {}
 LadderOpT(index_type const& i_idx, bool i_creatorF) : idx(i_idx), creatorF(i_creatorF)  {}
 LadderOpT& operator=(LadderOpT const& v)  { idx = v.idx;  creatorF = v.creatorF;  return *this; }
 
 bool operator< (LadderOpT const& v) const
  { if( creatorF != v.creatorF )  return( creatorF < v.creatorF );
    if( creatorF )  return idx_is_smaller( v.idx, idx );
    else            return idx_is_smaller( idx, v.idx ); }
 bool operator==(LadderOpT const& v) const  { return( creatorF == v.creatorF && idx_is_equal( idx, v.idx ) ); }
 bool operator!=(LadderOpT const& v) const  { return !( *this == v ); }
 bool operator> (LadderOpT const& v) const  { return !( *this < v || *this == v ); }
 
 index_type const& index() const  { return idx; }
 index_type&       index()        { return idx; }
 bool   is_valid()   const  { return idx != idx_invalid<index_type>(); }
 bool   is_creator() const  { return creatorF; }
 void   herm()	            { creatorF^=1; }
 
 PStream::opstream& write(PStream::opstream &out) const { out << idx << creatorF;  return out; }
 PStream::ipstream& read (PStream::ipstream &in)        { in  >> idx >> creatorF;  return in; }
 
 protected:
 index_type idx;
 bool       creatorF;
 };

inline LadderOp   ladder_p(sdp_gs::Index idx)  { return LadderOp(idx,1); }
inline LadderOp   ladder_m(sdp_gs::Index idx)  { return LadderOp(idx,0); }
inline LadderOpQn ladder_p(QuantumNumber const& qn)  { return LadderOpQn(qn,1); }
inline LadderOpQn ladder_m(QuantumNumber const& qn)  { return LadderOpQn(qn,0); }


//-----------------------------------------------------------------------[ stream stuff ]-
template <typename Idx>
inline PStream::opstream& operator<<(PStream::opstream& out, LadderOpT<Idx> const& v)
 { return v.write(out); }

template <typename Idx>
inline PStream::ipstream& operator>>(PStream::ipstream& in , LadderOpT<Idx>      & v)
 { return v.read (in);  }

template <typename Idx>
inline std::ostream& operator<<(std::ostream& out, LadderOpT<Idx> const& v)
 { out<<(v.is_creator()?"p":"m")<<v.index();  return out; }


//----------------------------------------------------------------------------[ algebra ]-
template <typename Idx>
inline void swap(LadderOpT<Idx>& v, LadderOpT<Idx>& w)  { LadderOpT<Idx> x(w);  w = v;  v = x; }

template <typename Idx>
inline LadderOpT<Idx> herm(LadderOpT<Idx> const& v)     { LadderOpT<Idx> ret(v);  ret.herm();  return ret; }


//---------------------------------------------------------------------------[ Distance ]-
template <typename T>
inline double Distance(LadderOpT<T> const& x, LadderOpT<T> const& y)
 { return std::max( Distance( x.index(), y.index() ),
                    Distance( x.is_creator(), y.is_creator() ) ); }



//__________________________________________________________________________________________________
//####################################################################################### Monomial #
template <typename T>
class MonomialT;

typedef MonomialT<LadderOp>   Monomial;
typedef MonomialT<LadderOpQn> MonomialQn;

typedef LinearAlgebra::Vector<Monomial>   Monomials;
typedef LinearAlgebra::Vector<MonomialQn> MonomialsQn;



//__________________________________________________________________________________________________
//####################################################################################### Monomial #
// An empty Monomial corresponds to value 1.
// CAUTION: A normal-ordered Monomial does not necessarily correspond to a normal-ordered MonomialQn (up to prefactor).
template <typename T>
class MonomialT
 {
 public:
 typedef LinearAlgebra::Vector<T> container_type;
 
 MonomialT()			     {}					// =1
 MonomialT(T const& v)		     { Op = container_type(1,v); }
 MonomialT(container_type const& v)  { Op = v; }
 MonomialT& operator=(MonomialT const& v)  { Op = v.Op;  return *this; }
 
 bool operator< (MonomialT const& v) const
  { if( size() != v.size() ) return( size() < v.size() );
    for(size_t n = 0; n < size(); n++)
     { if( Op[n] != v.Op[n] ) return( Op[n] < v.Op[n] ); }  return false; }
 bool operator==(MonomialT const& v) const  { return( Op == v.Op ); }
 bool operator!=(MonomialT const& v) const  { return !( *this == v ); }
 bool operator> (MonomialT const& v) const  { return !( *this < v || *this == v ); }
 
 size_t size() const  { return Op.size(); }
 
 container_type const& operator()() const  { return Op; }
 container_type      & operator()()        { return Op; }
 T const& operator[](size_t n) const  { return Op[n]; }
 T      & operator[](size_t n)        { return Op[n]; }
 
 MonomialT& operator*=(MonomialT const& v)  { Op = LinearAlgebra::direct_sum( Op, v.Op );  return *this; }
 bool is_valid() const  { for(size_t n=0; n<size(); n++) { if(!Op[n].is_valid()) return 0; }  return 1; }
 void herm()            { LinearAlgebra::mirror(Op);  for(size_t n=0; n<size(); n++) Op[n].herm(); }
 
 PStream::opstream& write(PStream::opstream &out) const  { out << Op;  return out; }
 PStream::ipstream& read (PStream::ipstream &in)         { in  >> Op;  return in; }
 
 protected:
 container_type Op;
 };


//-----------------------------------------------------------------------[ stream stuff ]-
template <typename T>
inline PStream::opstream& operator<<(PStream::opstream& out, MonomialT<T> const& v)
 { return v.write(out); }

template <typename T>
inline PStream::ipstream& operator>>(PStream::ipstream& in , MonomialT<T>      & v)
 { return v.read (in);  }

template <typename T>
std::ostream& operator<<(std::ostream& out, MonomialT<T> const& v);


//----------------------------------------------------------------------------[ algebra ]-
template <typename T>
MonomialT<T> operator*(MonomialT<T> const& v, MonomialT<T> const& w)
 { MonomialT<T> ret(v);  ret*=w;  return ret; }
inline Monomial   operator*(Monomial   const& v, Monomial   const& w)
 { Monomial ret(v);  ret*=w;  return ret; }
inline MonomialQn operator*(MonomialQn const& v, MonomialQn const& w)
 { MonomialQn ret(v);  ret*=w;  return ret; }
// template<>
// MonomialT<LadderOp> operator*<LadderOp>(MonomialT<LadderOp> const& v, MonomialT<LadderOp> const& w);
// template<>
// MonomialT<LadderOpQn> operator*<LadderOpQn>(MonomialT<LadderOpQn> const& v, MonomialT<LadderOpQn> const& w);

template <typename T>
inline MonomialT<T> herm(MonomialT<T> const& v)
 { MonomialT<T> ret(v);  ret.herm();  return ret; }

template <typename S>
inline Monomial
annihilatorProduct(S const& idcs, typename boost::enable_if<LinearAlgebra::is_vector<S> >::type* dummy = 0);


//----------------------------------------------------------------------[ quantumNumber ]-
QuantumNumber quantumNumber(MonomialQn const& x, sdp_gs::Indices const& qn_idcs);


//---------------------------------------------------------------------------[ Distance ]-
template <typename T>
inline double Distance(MonomialT<T> const& x, MonomialT<T> const& y);



//__________________________________________________________________________________________________
//##################################################################################### Polynomial #
template <typename T>
class PolynomialT;

typedef PolynomialT<LadderOp>   Polynomial;
typedef PolynomialT<LadderOpQn> PolynomialQn;

typedef LinearAlgebra::Vector<Polynomial>   Polynomials;
typedef LinearAlgebra::Vector<PolynomialQn> PolynomialsQn;



//__________________________________________________________________________________________________
//##################################################################################### Polynomial #
// CAUTION: A normal-ordered Polynomial does not necessarily correspond to a normal-ordered PolynomialQn (up to prefactor).
template <typename T>
class PolynomialT
 {
 public:
 typedef sdp_gs::CFloat coeff_type;
 typedef MonomialT<T>   monom_type;
 typedef std::map<monom_type,coeff_type>		container_type;
 typedef typename container_type::iterator		coeff_it;
 typedef typename container_type::const_iterator 	coeff_cit;
 typedef typename container_type::reverse_iterator	coeff_rit;
 
 PolynomialT()					{}				// =0
 PolynomialT(monom_type const& v)			{ Coeff[v] = 1; }
 PolynomialT(monom_type const& v, coeff_type c)	{ Coeff[v] = c; }		// =c*monom
 PolynomialT& operator=(PolynomialT const& v)	{ Coeff = v.Coeff;  return *this; }
 PolynomialT& operator=(monom_type const& v)	{ Coeff.clear();  Coeff[v] = 1;  return *this; }
 
 bool operator==(PolynomialT const& v)	{ return Coeff == v.Coeff; }
 
 size_t size() const				{ return Coeff.size(); }
 
 coeff_cit begin() const			{ return Coeff.begin(); }
 coeff_cit end()   const			{ return Coeff.end(); }
 coeff_it  begin()				{ return Coeff.begin(); }
 coeff_it  end()				{ return Coeff.end(); }
 coeff_rit rbegin()				{ return Coeff.rbegin(); }
 void      erase(coeff_it pos)			{ Coeff.erase(pos); }
 void      erase(monom_type const& v)		{ Coeff.erase(v); }
 
 coeff_type  operator[](monom_type const& v) const
  { coeff_cit i = Coeff.find(v);
    if( i == Coeff.end() ) return(0);  else return(i->second); }
 coeff_type& operator[](monom_type const& v)	{ return Coeff[v]; }
 
 PolynomialT& operator+=(PolynomialT const& v);
 PolynomialT& add(coeff_type c, PolynomialT const& v);	// optimization for the expression: +=c*v.
 PolynomialT& operator-=(PolynomialT const& v);
 PolynomialT  operator* (PolynomialT const& v) const;	//TODO use of erase_zero could be optimized
 PolynomialT& operator*=(coeff_type c);
 void erase_zeros()
  { coeff_it i = begin();
    while( i != end() )
     { if( i->second==coeff_type(0) ) erase(i++);  else ++i;  } }
 
 PStream::opstream& write(PStream::opstream &out) const  { out << Coeff;  return out; }
 PStream::ipstream& read (PStream::ipstream &in)         { in  >> Coeff;  return in; }

 protected:
 container_type Coeff;
 };


//-----------------------------------------------------------------------[ stream stuff ]-
template <typename T>
inline PStream::opstream& operator<<(PStream::opstream& out, PolynomialT<T> const& v)
 { return v.write(out); }

template <typename T>
inline PStream::ipstream& operator>>(PStream::ipstream& in , PolynomialT<T>      & v)
 { return v.read (in);  }

template <typename T>
std::ostream& operator<<(std::ostream& out, PolynomialT<T> const& v);


//----------------------------------------------------------------------------[ algebra ]-
inline PolynomialQn operator*(PolynomialQn const& v, MonomialQn const& w)  { return v*MonomialQn(w); }
inline PolynomialQn operator*(MonomialQn const& v, PolynomialQn const& w)  { return MonomialQn(v)*w; }

template <typename T>
inline PolynomialT<T> operator+(PolynomialT<T> const& v, PolynomialT<T> const& w)
 { PolynomialT<T> ret(v);  ret += w;  return ret; }
inline PolynomialQn operator+(PolynomialQn const& v, PolynomialQn const& w)
 { PolynomialQn ret(v);  ret += w;  return ret; }

template <typename T>
inline PolynomialT<T> operator-(PolynomialT<T> const& v, PolynomialT<T> const& w)
 { PolynomialT<T> ret(v);  ret -= w;  return ret; }
inline PolynomialQn operator-(PolynomialQn const& v, PolynomialQn const& w)
 { PolynomialQn ret(v);  ret -= w;  return ret; }

template <typename S, typename T>
inline typename boost::enable_if< LinearAlgebra::is_Fundamental<S>, PolynomialT<T> >::type
operator*(S const& c, MonomialT<T> const& v)
 { return PolynomialT<T>(v,c); }
template <typename S>
inline typename boost::enable_if< LinearAlgebra::is_Fundamental<S>, PolynomialQn >::type
operator*(S const& c, MonomialQn const& v)  { return PolynomialQn(v,c); }

template <typename S, typename T>
inline typename boost::enable_if< LinearAlgebra::is_Fundamental<S>, PolynomialT<T> >::type
operator*(S const& c, PolynomialT<T> const& v);

template <typename T>
PolynomialT<T> herm(PolynomialT<T> const& v);

template <typename T>
LinearAlgebra::Vector<PolynomialT<T> > herm(LinearAlgebra::Vector<PolynomialT<T> > const& v);


//----------------------------------------------------------------------[ quantumNumber ]-
QuantumNumber quantumNumber(PolynomialQn const& x, sdp_gs::Indices const& qn_idcs);

//-----------------------------------------------------------------------------[ degree ]-
template <typename T>
int degree(PolynomialT<T> const& v);

//-------------------------------------------------------------------------[ degree_min ]-
template <typename T>
int degree_min(PolynomialT<T> const& v);



//__________________________________________________________________________________________________
//######################################################################################### random #
inline LadderOp random_LadderOp(size_t NoModes)  { return LadderOp( rand()%NoModes, rand()%2 ); }

Monomial   random_Monomial  (size_t NoModes, size_t degree);

Polynomial random_Polynomial(size_t NoModes, sdp_gs::Indices const& NoTerms);


inline LadderOp random_LadderOp(sdp_gs::Indices const& subsystem)
 { return LadderOp( subsystem[ rand()%subsystem.size() ], rand()%2 ); }

Monomial random_Monomial(sdp_gs::Indices const& subsystem, size_t degree);

Polynomial random_Polynomial(LinearAlgebra::Vector<sdp_gs::Indices> const& subsystems,
                             sdp_gs::Indices const& NoTerms);



//##################################################################################################


#include "operator00.cc"

#endif // SDP_GS_KIT_OPERATOR00_H
