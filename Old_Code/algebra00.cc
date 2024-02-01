//##################################################################################################

//----------------------------------------------------------------------------[ is_zero ]-
template <typename T>
bool is_zero(MonomialT<T> const& v, int algebraF)
 {
 DEBUG_RANGE_CHECK( algebraF, 0, 2 );
 if( algebraF == 0 || !v.size() )  return false;
 
 for(size_t n = 0; n < v.size()-1; n++)
  for(size_t m = n+1; m < v.size(); m++)
   {
   if( v[m] == v[n] )  return true;
   if( v[m] == herm(v[n]) )  break;
   }
 
 return false;
 }


template <typename T>
bool is_zero(PolynomialT<T> const& v, int algebraF)
 {
 DEBUG_RANGE_CHECK( algebraF, 0, 2 );
 if( algebraF == 0 )  return false;
 
 for(typename PolynomialT<T>::coeff_cit i = v.begin(); i != v.end(); ++i)
  if( !( i->second == PolynomialT<T>::coeff_type(0) || is_zero( i->first, algebraF ) ) )  return false;
 
 return true;
 }


//------------------------------------------------------------------------[ erase_zeros ]-
template <typename T>
void erase_zeros(PolynomialT<T>& v, int algebraF)
 {
 typename PolynomialT<T>::coeff_it i = v.begin();
 while( i != v.end() )
  {
  if( i->second == PolynomialT<T>::coeff_type(0) || is_zero( i->first, algebraF ) )  v.erase(i++);
  else ++i;
  }
 }


//----------------------------------------------------------------------------[ commute ]-
template <typename T>
std::pair<int,int> commute(LadderOpT<T> const& op1, LadderOpT<T> const& op2, int algebraF)
 {
 if( algebraF == 2 )
  {
  if( idx_is_equal( op1.index(), op2.index() ) )  algebraF = 1;
  else  algebraF = 0;
  }
 switch( algebraF )
  {
  case 0:  if( op1 == herm(op2) )
            {
            if( op1.is_creator() )  return( std::pair<int,int>( 1, -1 ) );
            else                    return( std::pair<int,int>( 1,  1 ) );
            }
           return( std::pair<int,int>( 1,  0 ) );
  case 1:  if( op1 == op2 )        return( std::pair<int,int>( 0, 0 )   );
           if( op1 == herm(op2) )  return( std::pair<int,int>( -1,  1 ) );
           return( std::pair<int,int>( -1,  0 ) );
  default: PANIC("Unkown algebra type.")(algebraF);
  }
 return std::pair<int,int>();
 }


//-------------------------------------------------------------------[ normalPreordered ]-
template <typename T>
PolynomialT<T> normalPreordered(MonomialT<T> const& v, int algebraF)
 {
 using namespace LinearAlgebra;
 if( !v.size() )  return PolynomialT<T>(v);
 MonomialT<T>   w(v);
 PolynomialT<T> ret;
 int coeff = 1;
 bool anyChange;
 do
  {
  anyChange = 0;
  for(size_t n = 0; n < w.size()-1; n++)
   if( w[n] > w[n+1] )
    {
    std::pair<int,int> comm = commute( w[n], w[n+1], algebraF );
    if( comm.second != 0 )
     {
     MonomialT<T> reduced( direct_sum( w()[range(0,n)], w()[range(n+2,w.size())] ) );
     if( !is_zero( reduced, algebraF ) )
       ret += coeff*comm.second * reduced;
     }
    if( !comm.first )  return ret;
    coeff *= comm.first;
    swap( w()[n], w()[n+1] );
    anyChange = true;
    }
  if( is_zero( w, algebraF ) )  return ret;
  }
 while( anyChange );
 return( ret + coeff*w );
 }


//----------------------------------------------------------------------[ normalOrdered ]-
template <typename T>
PolynomialT<T> normalOrdered(PolynomialT<T> const& v, int algebraF)
 {
 RANGE_CHECK( algebraF, 0, 2 );
 PolynomialT<T> ret( v );
 for(int d = degree(v); d >= 2; d--)
  {
  PolynomialT<T> w;
  for(typename PolynomialT<T>::coeff_cit i = ret.begin(); i != ret.end(); ++i)
   {
   if( i->first.size() != d )  w += i->second * i->first;
   else  w.add( i->second, normalPreordered( i->first, algebraF ) );
   }
  ret = w;
  }
 return ret;
 }


//-----------------------------------------------------------------------[ is_hermitian ]-
template <typename T>
bool is_hermitian(PolynomialT<T> const& v, int algebraF)
 {
 PolynomialT<T> p = v - herm(v);
 return( normalOrdered(p,algebraF).size() == 0 );
 }


//-------------------------------------------------------------------[ is_normalOrdered ]-
template <typename T>
bool is_normalOrdered(MonomialT<T> const& v, int algebraF)
 {
 PolynomialT<T> p = normalOrdered( PolynomialT<T>(v), algebraF );
 if( p.size() != 1 )  return false;
 if( p.begin()->second != sdp_gs::CFloat(1) )  return false;
 return true;
 }


template <typename T>
bool is_normalOrdered(PolynomialT<T> const& v, int algebraF)
 {
 PolynomialT<T> p = normalOrdered( v, algebraF );
 return p == v;
 }


//##################################################################################################
