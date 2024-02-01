//__________________________________________________________________________________________________
//########################################################################################## idx_* #

//------------------------------------------------------------------------[ idx_invalid ]-
template <>
inline sdp_gs::Index  idx_invalid<sdp_gs::Index >()  { return std::numeric_limits<sdp_gs::Index>::max(); }

template <>
inline sdp_gs::PointF idx_invalid<sdp_gs::PointF>()  { return sdp_gs::PointF(); }


//-----------------------------------------------------------------------[ idx_is_equal ]-
template <>
inline bool idx_is_equal<sdp_gs::Index>(sdp_gs::Index const& idx1, sdp_gs::Index const& idx2)
 { return( idx1 == idx2 ); }

template <>
inline bool idx_is_equal<sdp_gs::PointF>(sdp_gs::PointF const& idx1, sdp_gs::PointF const& idx2)
 {
 if( idx1.size() != idx2.size() )  return false;
 LinearAlgebra::const_iterator<sdp_gs::PointF>::type idx1I = iterate(idx1);
 LinearAlgebra::const_iterator<sdp_gs::PointF>::type idx2I = iterate(idx2);
 for(; idx1I; ++idx1I, ++idx2I)
  if( Distance( *idx1I, *idx2I ) > 100*zeroThresh )  return false;
 return true;
 }


//---------------------------------------------------------------------[ idx_is_smaller ]-
template <>
inline bool idx_is_smaller<sdp_gs::Index>(sdp_gs::Index const& idx1, sdp_gs::Index const& idx2)
 { return( idx1 < idx2 ); }

template <>
inline bool idx_is_smaller<sdp_gs::PointF>(sdp_gs::PointF const& idx1, sdp_gs::PointF const& idx2)
 {
 if( idx1.size() != idx2.size() )  return( idx1.size() < idx2.size() );
 LinearAlgebra::const_iterator<sdp_gs::PointF>::type idx1I = iterate(idx1);
 LinearAlgebra::const_iterator<sdp_gs::PointF>::type idx2I = iterate(idx2);
 for(; idx1I; ++idx1I, ++idx2I)
  if( Distance( *idx1I, *idx2I ) > 100*zeroThresh )  return( *idx1I < *idx2I );
 return false;
 }


//##################################################################################################
