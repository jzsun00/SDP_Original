//__________________________________________________________________________________________________
//####################################################################################### Monomial #

//-----------------------------------------------------------------------[ stream stuff ]-
template <typename T>
std::ostream& operator<<(std::ostream& out, MonomialT<T> const& v)
 {
 if( v.size() ) out << v[0];
 else out << "Id";
 for(size_t n = 1; n < v.size(); n++)
  out << "*" << v[n];
 return out;
 }


//----------------------------------------------------------------------------[ algebra ]-
template <typename S>
inline Monomial
annihilatorProduct(S const& idcs, typename boost::enable_if<LinearAlgebra::is_vector<S> >::type* dummy)
 {
 Monomial::container_type ops( idcs.size() );
 for(size_t n = 0; n < idcs.size(); n++)  ops[n] = LadderOp( idcs[n], 0 );
 return Monomial( ops );
 }


//---------------------------------------------------------------------------[ Distance ]-
template <typename T>
inline double Distance(MonomialT<T> const& x, MonomialT<T> const& y)
 {
 if( x.size() != y.size() )  return Distance( x.size(), y.size() );
 double ret = 0;
 for(size_t n = 0; n < x.size(); n++)  ret = std::max( ret, Distance( x[n], y[n] ) );
 return ret;
 }



//__________________________________________________________________________________________________
//##################################################################################### Polynomial #

//-----------------------------------------------------------------------[ stream stuff ]-
template <typename T>
std::ostream& operator<<(std::ostream& out, PolynomialT<T> const& v)
 {
 typename PolynomialT<T>::coeff_cit i = v.begin();
 
 if( i != v.end() )
  {
  if( sdp_gs::imag(i->second) == 0 )  out << sdp_gs::real(i->second) << "*" << i->first;
  else  out << i->second << "*" << i->first;
  ++i;
  }
 else out << "0";
 
 for(; i != v.end(); ++i )
  {
  if( sdp_gs::imag(i->second) == 0 )
   out << (sdp_gs::real(i->second) >= 0? " + ":" - " ) << std::abs(sdp_gs::real(i->second)) << "*" << i->first;
  else
   out << " + " << i->second << "*" << i->first;
  }
 
 return out;
 }


//-------------------------------------------------------------------------[ operator+= ]-
// REMARK: One could do the sum slightly more efficiently for the case, where the "supports" of
//          the polynomials are comparable.
template <typename T>
PolynomialT<T>& PolynomialT<T>::operator+=(PolynomialT const& v)
 {
 for(coeff_cit i = v.begin(); i != v.end(); ++i)
  {
  coeff_it j = Coeff.find(i->first);
  if( j == Coeff.end() )  Coeff.insert(*i);
  else
   {
   if( j->second + i->second == coeff_type(0) )  Coeff.erase(j);
   else  j->second += i->second;
   }
  }
 return *this;
 }


//--------------------------------------------------------------------------------[ add ]-
template <typename T>
PolynomialT<T>& PolynomialT<T>::add(coeff_type c, PolynomialT const& v)
 {
 if( c == coeff_type(0) || !v.size() )  return *this;
 for(coeff_cit i = v.begin(); i != v.end(); ++i)
  {
  coeff_it j = Coeff.find(i->first);
  if( j == Coeff.end() )  Coeff[i->first] = c*i->second;
  else
   {
   if( j->second + c*i->second == coeff_type(0) )  Coeff.erase(j);
   else  j->second += c*i->second;
   }
  }
 return *this;
 }


//-------------------------------------------------------------------------[ operator-= ]-
template <typename T>
PolynomialT<T>& PolynomialT<T>::operator-=(PolynomialT const& v)
 {
 for(coeff_cit i = v.begin(); i != v.end(); ++i)
  {
  coeff_it j = Coeff.find(i->first);
  if( j == Coeff.end() )  Coeff[i->first] = -i->second;
  else
   {
   if( j->second - i->second == coeff_type(0) )  Coeff.erase(j);
   else  j->second -= i->second;
   }
  }
 return *this;
 }


//--------------------------------------------------------------------------[ operator* ]-
template <typename T>
PolynomialT<T> PolynomialT<T>::operator*(PolynomialT const& v) const
 {
 PolynomialT ret;
 for(coeff_cit i = begin(); i != end(); ++i)
  for(coeff_cit j = v.begin(); j != v.end(); ++j)
   ret[i->first*j->first] += i->second*j->second;
 ret.erase_zeros();
 return ret;
 }


//-------------------------------------------------------------------------[ operator*= ]-
template <typename T>
PolynomialT<T>& PolynomialT<T>::operator*=(coeff_type c)
 {
 if( c == coeff_type(0) )  (*this) = PolynomialT();
 else
  for(coeff_it i = begin(); i != end(); ++i)
   i->second *= c;
 return *this;
 }


//--------------------------------------------------------------------------[ operator* ]-
template <typename S, typename T>
inline typename boost::enable_if< LinearAlgebra::is_Fundamental<S>, PolynomialT<T> >::type
operator*(S const& c, PolynomialT<T> const& v)
 {
 if(    (typename PolynomialT<T>::coeff_type)c
     == (typename PolynomialT<T>::coeff_type)0 )  return PolynomialT<T>();
 PolynomialT<T> ret(v);
 ret *= c;
 return ret;
 }


//-------------------------------------------------------------------------------[ herm ]-
//assuming that "herm" is injective in the set of all monomes
template <typename T>
PolynomialT<T> herm(PolynomialT<T> const& v)
 {
 PolynomialT<T> ret;
 for(typename PolynomialT<T>::coeff_cit i = v.begin(); i != v.end(); ++i)
  ret[::herm(i->first)] = LinearAlgebra::conj(i->second);
 return ret;
 }


//-------------------------------------------------------------------------------[ herm ]-
template <typename T>
LinearAlgebra::Vector<PolynomialT<T> > herm(LinearAlgebra::Vector<PolynomialT<T> > const& v)
 {
 LinearAlgebra::Vector<PolynomialT<T> > ret( v.size() );
 for(size_t n = 0; n < v.size(); n++)  ret[n] = herm(v[n]);
 return ret;
 }


//-----------------------------------------------------------------------------[ degree ]-
template <typename T>
int degree(PolynomialT<T> const& v)
 {
 int r = -1;
 for(typename PolynomialT<T>::coeff_cit i = v.begin(); i != v.end(); ++i)
  r = std::max( r, (int)i->first.size() );
 return r;
 }

//-------------------------------------------------------------------------[ degree_min ]-
template <typename T>
int degree_min(PolynomialT<T> const& v)
 {
 if( !v.size() )  return -1;
 typename PolynomialT<T>::coeff_cit i = v.begin();
 int r = (int)i->first.size();
 for(++i; i != v.end(); ++i)
  r = std::min( r, (int)i->first.size() );
 return r;
 }



//##################################################################################################
