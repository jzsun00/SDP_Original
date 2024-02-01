//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 29.06.2011 Thomas Barthel
//% SDP groundstate kit
//% Algebra for operators in second quantization.
//% Supported algebraF values:
//%  0: bosons        a_i a_j^+ - a_j^+ a_i = \delta_{ij},    a_i a_j =  a_j a_i
//%  1: fermions      a_i a_j^+ + a_j^+ a_i = \delta_{ij},    a_i a_j = -a_j a_i
//%  2: hrdcr. bosons a_i a_j^+ - (-1)^{\delta_{ij}} a_j^+ a_i = \delta_{ij},  a_i a_j = \delta_{i!=j}a_j a_i
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// key words: TODO, QUESTION, TEST, CHECK, REMARK, CAUTION
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#ifndef SDP_GS_KIT_ALGEBRA00_H
#define SDP_GS_KIT_ALGEBRA00_H

#include "operator00.h"



//##################################################################################################

//----------------------------------------------------------------------------[ is_zero ]-
template <typename T>
bool is_zero(MonomialT<T>   const& v, int algebraF);

template <typename T>
bool is_zero(PolynomialT<T> const& v, int algebraF);


//------------------------------------------------------------------------[ erase_zeros ]-
template <typename T>
void erase_zeros(PolynomialT<T>& v, int algebraF);

template <typename T>
inline PolynomialT<T> erased_zeros(PolynomialT<T> const& v, int algebraF)
 { PolynomialT<T> ret(v);  erase_zeros( ret, algebraF );  return ret; }


//----------------------------------------------------------------------------[ commute ]-
// op1*op2 = s*op2*op1 + c;  The return value of the function is (s,c).
template <typename T>
std::pair<int,int> commute(LadderOpT<T> const& op1, LadderOpT<T> const& op2, int algebraF);


//-------------------------------------------------------------------[ normalPreordered ]-
template <typename T>
PolynomialT<T> normalPreordered(MonomialT<T> const& v, int algebraF);


//----------------------------------------------------------------------[ normalOrdered ]-
template <typename T>
PolynomialT<T> normalOrdered(PolynomialT<T> const& v, int algebraF);


//-----------------------------------------------------------------------[ is_hermitian ]-
template <typename T>
bool is_hermitian(PolynomialT<T> const& v, int algebraF);


//-------------------------------------------------------------------[ is_normalOrdered ]-
template <typename T>
bool is_normalOrdered(MonomialT<T> const& v, int algebraF);

template <typename T>
bool is_normalOrdered(PolynomialT<T> const& v, int algebraF);



//##################################################################################################


#include "algebra00.cc"

#endif // SDP_GS_KIT_ALGEBRA00_H
