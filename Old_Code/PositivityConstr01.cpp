#include "PositivityConstr01.h"



//__________________________________________________________________________________________________
//########################################################################## PositivityConstraints #

//-----------------------------------------------------------------------[ stream stuff ]-
template<>
const unsigned sdp_gs::ClassID<PositivityConstraintsF >::ID = 7553901;
template<>
const unsigned sdp_gs::ClassID<PositivityConstraintsCF>::ID = 7553902;


//---------------------------------------------------------------------[ from_rootBasis ]-
// Expanding products s^+*s' for all s, s' from an operator basis in the Green's function basis.
template<>
PositivityConstraintsCF
PositivityConstraintsCF::from_rootBasis(PolynomialsQn const& opBasis, GreenBasis const& G_basis, double ZeroThresh)
 {
 PositivityConstraintsCF ret;
 using namespace LinearAlgebra;
 ret.M_Re = data_type( G_basis.sizeReal(), matrix_type( opBasis.size(), opBasis.size() ) );
 ret.M_Im = data_type( G_basis.sizeImag(), matrix_type( opBasis.size(), opBasis.size() ) );
 TRACE( opBasis.size() );
 size_t m0 = 1;
 for(size_t m = 0; m < opBasis.size(); m++)
  {
  if( m == m0 )  { std::cout << " " << m << std::flush;  m0*=2; }
  for(size_t n = 0; n < opBasis.size(); n++)
   {
   GreenFunctionDual P = G_basis.expand_dual( herm(opBasis[m])*opBasis[n] );
   for(size_t s = 0; s < G_basis.sizeReal(); s++)
    if( std::abs( P.first[s] )  > ZeroThresh )  ret.M_Re[s](m,n) = P.first[s];
   for(size_t s = 0; s < G_basis.sizeImag(); s++)
    if( std::abs( P.second[s] ) > ZeroThresh )  ret.M_Im[s](m,n) = sdp_gs::my_i*P.second[s];
   }
  }
 #if !defined(NDEBUG)
 for(size_t s = 0; s < G_basis.sizeReal(); s++)
  { CHECK_COMPARE( normInf( ret.M_Re[s] - herm(ret.M_Re[s]) ), <, zeroThresh ); }
 for(size_t s = 0; s < G_basis.sizeImag(); s++)
  { CHECK_COMPARE( normInf( ret.M_Im[s] - herm(ret.M_Im[s]) ), <, zeroThresh ); }
 #endif
 return ret;
 }


//--------------------------------------------------------------------[ from_observable ]-
// Expanding an hermitian operator.
template<>
PositivityConstraintsF
PositivityConstraintsF::from_observable(PolynomialQn const& op, GreenBasis const& G_basis, double ZeroThresh)
 {
 PositivityConstraintsF ret;
 using namespace LinearAlgebra;
 ret.M_Re = data_type( G_basis.sizeReal(), matrix_type( 1, 1 ) );
 ret.M_Im = data_type( G_basis.sizeImag(), matrix_type( 1, 1 ) );
 GreenFunctionDualHerm P = greenFunctionDualHerm( G_basis.expand_dual( op ) );
 for(size_t s = 0; s < G_basis.sizeReal(); s++)
  if( std::abs( P.first[s] )  > ZeroThresh )  ret.M_Re[s](0,0) = P.first[s];
 for(size_t s = 0; s < G_basis.sizeImag(); s++)
  if( std::abs( P.second[s] ) > ZeroThresh )  ret.M_Im[s](0,0) = P.second[s];
 return ret;
 }


//-------------------------------------------------------------------------------[ real ]-
PositivityConstraintsF real(PositivityConstraintsCF const& constr)
 {
 using namespace LinearAlgebra;
 PositivityConstraintsF ret;
 size_t D = constr.dim();
 
 if( constr.nnz_Im_max() == 0 && constr.normInf_imag_Re() < 100*zeroThresh )
  return all_real( constr );
 
 ret.matricesReal() = PositivityConstraintsF::data_type( constr.sizeReal() );
 for(size_t n = 0; n < constr.sizeReal(); n++)
  {
  PositivityConstraintsF::matrix_type M( 2*D, 2*D );
  set_sparse_block( M, 0, 0,  real( constr.matricesReal()[n] ) );
  set_sparse_block( M, 0, D, -imag( constr.matricesReal()[n] ) );
  set_sparse_block( M, D, 0,  imag( constr.matricesReal()[n] ) );
  set_sparse_block( M, D, D,  real( constr.matricesReal()[n] ) );
//   M( range(0,  D), range(0,  D) ) =  real( constr.matricesReal()[n] );
//   M( range(0,  D), range(D,2*D) ) = -imag( constr.matricesReal()[n] );
//   M( range(D,2*D), range(0,  D) ) =  imag( constr.matricesReal()[n] );
//   M( range(D,2*D), range(D,2*D) ) =  real( constr.matricesReal()[n] );
  ret.matricesReal()[n] = M;
  }
 
 ret.matricesImag() = PositivityConstraintsF::data_type( constr.sizeImag() );
 for(size_t n = 0; n < constr.sizeImag(); n++)
  {
  PositivityConstraintsF::matrix_type M( 2*D, 2*D );
  set_sparse_block( M, 0, 0,  real( constr.matricesImag()[n] ) );
  set_sparse_block( M, 0, D,  imag( constr.matricesImag()[n] ) );
  set_sparse_block( M, D, 0, -imag( constr.matricesImag()[n] ) );
  set_sparse_block( M, D, D,  real( constr.matricesImag()[n] ) );
/*  set_sparse_block( M, 0, 0, -imag( constr.matricesImag()[n] ) );
  set_sparse_block( M, 0, D, -real( constr.matricesImag()[n] ) );
  set_sparse_block( M, D, 0,  real( constr.matricesImag()[n] ) );
  set_sparse_block( M, D, D, -imag( constr.matricesImag()[n] ) );*/
//   M( range(0,  D), range(0,  D) ) = -imag( constr.matricesImag()[n] );
//   M( range(0,  D), range(D,2*D) ) = -real( constr.matricesImag()[n] );
//   M( range(D,2*D), range(0,  D) ) =  real( constr.matricesImag()[n] );
//   M( range(D,2*D), range(D,2*D) ) = -imag( constr.matricesImag()[n] );
  ret.matricesImag()[n] = M;
  }
 
 return ret;
 }


//---------------------------------------------------------------------------[ all_real ]-
PositivityConstraintsF all_real(PositivityConstraintsCF const& constr)
 {
 using namespace LinearAlgebra;
 PRECONDITION( constr.nnz_Im_max() == 0 );
 CHECK_COMPARE( constr.normInf_imag_Re(), <, 100*zeroThresh );
 PositivityConstraintsF ret;
 
 ret.matricesReal() = PositivityConstraintsF::data_type( constr.sizeReal() );
 for(size_t n = 0; n < constr.sizeReal(); n++)
  ret.matricesReal()[n] = real( constr.matricesReal()[n] );
 
 return ret;
 }



//##################################################################################################
