//##################################################################################################

//-----------------------------------------------------------------------[ distance_min ]-
template <typename T, typename L>
double distance_min(T const& x, L const& list,
  typename boost::enable_if< LinearAlgebra::is_vector<L> >::type* dummy1,
  typename boost::enable_if< boost::is_same<T, typename LinearAlgebra::interface<L>::value_type> >::type* dummy2)
 {
 double dist = 3e8;
 for(size_t n = 0; n < list.size(); n++)
  dist = std::min( dist, Distance( list[n], x ) );
 return dist;
 }

template <typename T, typename L>
double distance_min(T const& x, L const& list, size_t& closest_element,
  typename boost::enable_if< LinearAlgebra::is_vector<L> >::type* dummy1,
  typename boost::enable_if< boost::is_same<T, typename LinearAlgebra::interface<L>::value_type> >::type* dummy2)
 {
 double dist = 3e8;
 closest_element = std::numeric_limits<size_t>::max();
 for(size_t n = 0; n < list.size(); n++)
  {
  double d = Distance( list[n], x );
  if( d < dist )
   {
   dist = d;
   closest_element = n;
   }
  }
 return dist;
 }

template <typename L1, typename L2>
double distance_min(L1 const& list1, L2 const& list2,
  typename boost::enable_if< LinearAlgebra::is_vector<L1> >::type* dummy1,
  typename boost::enable_if< LinearAlgebra::is_vector<L2> >::type* dummy2,
  typename boost::enable_if< boost::is_same<typename LinearAlgebra::interface<L1>::value_type,
                                            typename LinearAlgebra::interface<L2>::value_type> >::type* dummy3)
 {
 double dist = 3e8;
 for(size_t n = 0; n < list1.size(); n++)
  dist = std::min( dist, distance_min( list1[n], list2 ) );
 return dist;
 }

template <typename L>
double distance_min(L const& list,
  typename boost::enable_if< LinearAlgebra::is_vector<L> >::type* dummy1)
 {
 double dist = 3e8;
 for(size_t n = 0; n < list.size(); n++)
  for(size_t m = n+1; m < list.size(); m++)
   dist = std::min( dist, Distance( list[n], list[m] ) );
 return dist;
 }



//__________________________________________________________________________________________________
//####################################################################################### GroupRep #
template <typename element_type>
template <typename T>
GroupRep<element_type>::
GroupRep(T const& i_elements, typename boost::enable_if<LinearAlgebra::is_vector<T> >::type* dummy)
 : Elements(i_elements)
 {
 for(size_t n = 0; n < Elements.size(); n++)
  {
  CHECK_EQUAL( Elements[n].size1(), Elements[n].size2() );
  CHECK_EQUAL( Elements[n].size1(), Elements[0].size1() );
  }
 }


//---------------------------------------------------------------------[ generate_group ]-
template <typename element_type>
template <typename T>
inline void GroupRep<element_type>::
generate(T const& generators, element_type const& Id,
         typename boost::enable_if<LinearAlgebra::is_vector<T> >::type* dummy)
 {
 using namespace LinearAlgebra;
 Elements = container_type( 1, Id );
 if( !generators.size() )  return;
 
 size_t N0 = 0;
 size_t N;
 do
  {
  N = Elements.size();
  for(size_t m = 0; m < generators.size(); m++)
   for(size_t n = N0; n < N; n++)
    if( distance_min( eval(Elements[n]*generators[m]), Elements ) > 100*zeroThresh )
     Elements << Elements[n]*generators[m];
  N0 = N;
  }
 while( Elements.size() > N );
 }



//__________________________________________________________________________________________________
//######################################################################### common symmetry groups #

//---------------------------------------------------------------------[ generate_group ]-
inline GroupRep<sdp_gs::MatrixF> generate_group(LinearAlgebra::Vector<sdp_gs::MatrixF> const& generators,
                                                sdp_gs::MatrixF const& Id)
 {
 GroupRep<sdp_gs::MatrixF> ret;
 ret.generate( generators, Id );
 return ret;
 }



//##################################################################################################
