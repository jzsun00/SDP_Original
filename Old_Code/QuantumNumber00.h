//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 04.07.2011 Thomas Barthel
//% SDP groundstate kit
//% Mode quantum numbers.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// key words: TODO, QUESTION, TEST, CHECK, REMARK, CAUTION
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#ifndef SDP_GS_KIT_QUANTUMNUMBER00_H
#define SDP_GS_KIT_QUANTUMNUMBER00_H

#include "configuration.h"

#include "linalg/LabelListL00.h"

#include "linearalgebra/pstreamio.h"



//__________________________________________________________________________________________________
//########################################################################################## idx_* #

//------------------------------------------------------------------------[ idx_invalid ]-
template <typename T>
inline T idx_invalid();

//-----------------------------------------------------------------------[ idx_is_equal ]-
template <typename T>
inline bool idx_is_equal(T const& idx1, T const& idx2);

//---------------------------------------------------------------------[ idx_is_smaller ]-
template <typename T>
inline bool idx_is_smaller(T const& idx1, T const& idx2);

//---------------------------------------------------------------------------[ Idx_less ]-
template <typename T>
struct Idx_less : std::binary_function<T,T,bool>
 { bool operator()(T const& idx1, T const& idx2) const  { return idx_is_smaller(idx1,idx2); } };



//__________________________________________________________________________________________________
//################################################################################## QuantumNumber #
typedef LinearAlgebra::Vector<sdp_gs::Float> QuantumNumber;
typedef LinearAlgebra::Vector<QuantumNumber> QuantumNumbers;
typedef LabelList<QuantumNumber,Idx_less<QuantumNumber> > OneParticleBasis;

typedef sdp_gs::Indices Subsystem;


//-----------------------------------------------------------------------[ qn_projected ]-
QuantumNumbers qn_projected(QuantumNumbers const& qns, sdp_gs::Indices const& idcs);

//-----------------------------------------------------------------------[ qn_directSum ]-
QuantumNumbers qn_directSum(QuantumNumbers const& qns1, QuantumNumbers const& qns2);




//##################################################################################################


#include "QuantumNumber00.cc"

#endif // SDP_GS_KIT_QUANTUMNUMBER00_H
