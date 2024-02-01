//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 22.07.2011 Thomas Barthel
//% SDP groundstate kit
//% Positivity constraints.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// key words: TODO, QUESTION, TEST, CHECK, REMARK, CAUTION
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#ifndef SDP_GS_KIT_POSITIVITYCONSTR01_H
#define SDP_GS_KIT_POSITIVITYCONSTR01_H

#include "algebra/VerySparseMatrix00.h"
#include "latticeAlgebra/GreenFunction00.h"



//##################################################################################################

//-------------------------------------------------------------------[ set_sparse_block ]-
template <typename T>
void set_sparse_block(LinearAlgebra::VerySparseMatrix<T>& M, sdp_gs::Index row0, sdp_gs::Index col0,
                      LinearAlgebra::VerySparseMatrix<T> const& m);


//-------------------------------------------------------------------[ get_sparse_block ]-
template <typename T>
LinearAlgebra::VerySparseMatrix<T>
get_sparse_block(LinearAlgebra::VerySparseMatrix<T> const& M, sdp_gs::Indices const& rows, sdp_gs::Indices const& cols);


//--------------------------------------------------------------------------------[ nnz ]-
template <typename T>
size_t nnz(LinearAlgebra::VerySparseMatrix<T> const& M)  { return M.size(); }



//__________________________________________________________________________________________________
//########################################################################## PositivityConstraints #
// Positivity constraints, resulting from G[C^+ C] >= 0 for all C from some operator vector space.
// With some basis {s} for this vector space, the operators C can be expanded in the form
//  C = \sum_s C_s * s and the constraints attain the form
//  G[C^+ C] = \sum_{s's''} C^*_{s'} G[(s')^+ s''] C_{s''}
//           = \sum_s\sum_{s's''} G_s C^*_{s'} [M_s]_{s',s''} C_{s''}  for all C
// <=> \sum_s G_s M_s = \sum_s (Re G_s + i Im G_s) M_s =  >= 0  (">=" meaning positive definite)
template <typename T>
class PositivityConstraints
 {
 public:
 typedef T value_type;
 typedef LinearAlgebra::VerySparseMatrix<T> matrix_type;
 typedef LinearAlgebra::Vector<matrix_type> data_type;
 
 PositivityConstraints()  {}
 PositivityConstraints(size_t N_Re, size_t N_Im) : M_Re(N_Re), M_Im(N_Im)  {}
 PositivityConstraints& operator=(PositivityConstraints const& v)
  { M_Re = v.M_Re;  M_Im = v.M_Im;  return *this; }
 
 size_t sizeReal() const  { return M_Re.size(); }
 size_t sizeImag() const  { return M_Im.size(); }
 size_t size()     const  { return( sizeReal() + sizeImag() ); }
 size_t dim() const;
 
 data_type const& matricesReal() const  { return M_Re; }
 data_type const& matricesImag() const  { return M_Im; }
 data_type&       matricesReal()        { return M_Re; }
 data_type&       matricesImag()        { return M_Im; }
 
 sdp_gs::Indices nnz_Re() const;
 sdp_gs::Indices nnz_Im() const;
 size_t          nnz_Im_max() const
  { sdp_gs::Indices nnz = nnz_Im();  return( nnz.size()==0? 0 : LinearAlgebra::max(nnz) ); }
 double normInf_imag_Re() const;
 double normInf_imag_Im() const;
 double normInf_real_Im() const;
 
 matrix_type eval(GreenFunction const& G) const;
 
 PStream::opstream& write(PStream::opstream &out) const;
 PStream::ipstream& read (PStream::ipstream &in);
 
 static PositivityConstraints<T>
 from_rootBasis(PolynomialsQn const& opBasis, GreenBasis const& G_basis,
                double ZeroThresh = zeroThresh);
 static PositivityConstraints<T>
 from_observable(PolynomialQn const& op, GreenBasis const& G_basis, double ZeroThresh = zeroThresh);
 
 protected:
 data_type M_Re;
 data_type M_Im;
 };


//---------------------------------------------------------------------------[ typedefs ]-
typedef PositivityConstraints<sdp_gs::Float>  PositivityConstraintsF;
typedef PositivityConstraints<sdp_gs::CFloat> PositivityConstraintsCF;


//-----------------------------------------------------------------------[ stream stuff ]-
template <typename T>
inline PStream::opstream& operator<<(PStream::opstream& out, PositivityConstraints<T> const& v)
 { return v.write(out); }
template <typename T>
inline PStream::ipstream& operator>>(PStream::ipstream& in , PositivityConstraints<T>      & v)
 { return v.read (in);  }
template <typename T>
std::ostream& operator<<(std::ostream& out, PositivityConstraints<T> const& v);



//-------------------------------------------------------------------------------[ real ]-
// Translate complex N x N constraints into real 2N x 2N constraints.
PositivityConstraintsF real(PositivityConstraintsCF const& constr);

inline PositivityConstraintsF real(PositivityConstraintsF const& constr)  { return constr; }


//---------------------------------------------------------------------------[ all_real ]-
// Translate complex N x N constraints into real N x N constraints under the assumption
//  that all constraint matrices are real and that the Hamiltonian is real -> real GF.
PositivityConstraintsF all_real(PositivityConstraintsCF const& constr);

inline PositivityConstraintsF all_real(PositivityConstraintsF const& constr)  { return constr; }


//---------------------------------------------------------------------[ blockStructure ]-
template <typename T>
LinearAlgebra::Vector<sdp_gs::Indices>
blockStructure(PositivityConstraints<T> const& constr);


//--------------------------------------------------------------------[ separate_blocks ]-
template <typename T>
LinearAlgebra::Vector<PositivityConstraints<T> >
separate_blocks(PositivityConstraints<T> const& constr, LinearAlgebra::Vector<sdp_gs::Indices> const& blockIdcs);


//----------------------------------------------------------------[ GreenBasisStructure ]-
// Determines groups of Green's basis elements that are coupled to the same set of constraint blocks.
// The function returns the list of groups and the corresponding sets of constraint indices.
template <typename T>
LinearAlgebra::Vector<sdp_gs::Indices>
GreenBasisStructure_Re(LinearAlgebra::Vector< PositivityConstraints<T> > const& constr,
                       LinearAlgebra::Vector<sdp_gs::Indices>& constraintIdcs);

template <typename T>
LinearAlgebra::Vector<sdp_gs::Indices>
GreenBasisStructure_Im(LinearAlgebra::Vector< PositivityConstraints<T> > const& constr,
                       LinearAlgebra::Vector<sdp_gs::Indices>& constraintIdcs);




//##################################################################################################


#include "PositivityConstr01.cc"

#endif // SDP_GS_KIT_POSITIVITYCONSTR01_H
