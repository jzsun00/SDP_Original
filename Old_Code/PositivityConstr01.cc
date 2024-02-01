//##################################################################################################

//-------------------------------------------------------------------[ set_sparse_block ]-
template <typename T>
void set_sparse_block(LinearAlgebra::VerySparseMatrix<T>& M, sdp_gs::Index row0, sdp_gs::Index col0,
                      LinearAlgebra::VerySparseMatrix<T> const& m)
 {
 using namespace LinearAlgebra;
 typedef LinearAlgebra::VerySparseMatrix<T> matrix_type;
 for(typename matrix_type::const_iterator i = m.begin(); i != m.end(); ++i)
   M( i->first.first + row0, i->first.second + col0 ) = i->second;
 }


//-------------------------------------------------------------------[ get_sparse_block ]-
template <typename T>
LinearAlgebra::VerySparseMatrix<T>
get_sparse_block(LinearAlgebra::VerySparseMatrix<T> const& M, sdp_gs::Indices const& rows, sdp_gs::Indices const& cols)
 {
 LabelList<sdp_gs::Index> rowIdx( rows );
 LabelList<sdp_gs::Index> colIdx( cols );
 
 using namespace LinearAlgebra;
 typedef VerySparseMatrix<T> matrix_type;
 matrix_type ret( rows.size(), cols.size() );
 
 for(typename matrix_type::const_iterator i = M.begin(); i != M.end(); ++i)
  if( rowIdx.contains( i->first.first ) && colIdx.contains( i->first.second ) )
   ret( rowIdx( i->first.first ), colIdx( i->first.second ) ) = i->second;
 
 return ret;
 }



//__________________________________________________________________________________________________
//########################################################################## PositivityConstraints #

//-----------------------------------------------------------------------[ stream stuff ]-
template <typename T>
std::ostream& operator<<(std::ostream& out, PositivityConstraints<T> const& v)
 {
 out << "(matricesReal=" << v.matricesReal() << ",\n"
     << " matricesImag=" << v.matricesImag()
     << ")";
 return out;
 }

template <typename T>
PStream::opstream& PositivityConstraints<T>::write(PStream::opstream &out) const
 {
 out << (unsigned) sdp_gs::ClassID<PositivityConstraints<T> >::ID;	// ID
 out << (int)1; 							// version
 out << M_Re << M_Im;
 return out;
 }

template <typename T>
PStream::ipstream& PositivityConstraints<T>::read(PStream::ipstream &in)
 {
 unsigned ID;
 in >> ID;
 CHECK_EQUAL( ID, sdp_gs::ClassID<PositivityConstraints<T> >::ID );
 int ver;
 in >> ver;
 DEBUG_CHECK_EQUAL(ver,1);
 in >> M_Re >> M_Im;
 return in;
 }


//--------------------------------------------------------------------------------[ dim ]-
template <typename T>
size_t PositivityConstraints<T>::dim() const
 {
 if( sizeReal() )  { DEBUG_CHECK_EQUAL( M_Re[0].size1(), M_Re[0].size2() );  return M_Re[0].size1(); }
 if( sizeImag() )  { DEBUG_CHECK_EQUAL( M_Im[0].size1(), M_Im[0].size2() );  return M_Im[0].size1(); }
 return 0;
 }


//-------------------------------------------------------------------------------[ eval ]-
template <typename T>
typename PositivityConstraints<T>::matrix_type
PositivityConstraints<T>::eval(GreenFunction const& G) const
 {
 matrix_type ret( dim(), dim() );
 for(size_t n = 0; n < sizeReal(); n++)
  ret += G.first[n]*M_Re[n];
 for(size_t n = 0; n < sizeImag(); n++)
  ret += G.second[n]*M_Im[n];
 return ret;
 }


//--------------------------------------------------------------------------------[ nnz ]-
template <typename T>
sdp_gs::Indices PositivityConstraints<T>::nnz_Re() const
 {
 sdp_gs::Indices ret( sizeReal() );
 for(size_t n = 0; n < sizeReal(); n++)
  ret[n] = nnz( matricesReal()[n] );
 return ret;
 }

template <typename T>
sdp_gs::Indices PositivityConstraints<T>::nnz_Im() const
 {
 sdp_gs::Indices ret( sizeImag() );
 for(size_t n = 0; n < sizeImag(); n++)
  ret[n] = nnz( matricesImag()[n] );
 return ret;
 }


//-----------------------------------------------------------------------[ normInf_imag ]-
template <typename T>
double PositivityConstraints<T>::normInf_imag_Re() const
 {
 using namespace LinearAlgebra;
 double ret = 0;
 for(size_t n = 0; n < sizeReal(); n++)
  ret = std::max( ret, normInf( imag(matricesReal()[n]) ) );
 return ret;
 }

template <typename T>
double PositivityConstraints<T>::normInf_imag_Im() const
 {
 using namespace LinearAlgebra;
 double ret = 0;
 for(size_t n = 0; n < sizeImag(); n++)
  ret = std::max( ret, normInf( imag(matricesImag()[n]) ) );
 return ret;
 }


//-----------------------------------------------------------------------[ normInf_real ]-
template <typename T>
double PositivityConstraints<T>::normInf_real_Im() const
 {
 using namespace LinearAlgebra;
 double ret = 0;
 for(size_t n = 0; n < sizeImag(); n++)
  ret = std::max( ret, normInf( real(matricesImag()[n]) ) );
 return ret;
 }


//---------------------------------------------------------------------[ blockStructure ]-
// If necessary, the following identification of the block structure could be done more efficiently.
namespace Private_PosConstr
{
template <typename T>
void update_blockIDs(sdp_gs::Indices& blockID, sdp_gs::Index& maxBlock,
                     typename PositivityConstraints<T>::matrix_type const& M)
 {
 using namespace LinearAlgebra;
 typedef typename PositivityConstraints<T>::matrix_type matrix_type;
 for(typename matrix_type::const_iterator i = M.begin(); i != M.end(); ++i)
  {
  sdp_gs::Index idx1 = i->first.first;
  sdp_gs::Index idx2 = i->first.second;
  if( blockID[idx1] > blockID[idx2] )  std::swap( idx1, idx2 );
  if( !blockID[idx1] )
   {
   if( !blockID[idx2] )		// both in no block => new block
    {
    maxBlock++;
    blockID[idx1] = blockID[idx2] = maxBlock;
    }
   else				// idx2 is in a block, idx1 is not => put idx1 in block[idx2]
    blockID[idx1] = blockID[idx2];
   }
  else if( blockID[idx1] != blockID[idx2] )	// both in separate blocks => unite the blocks
   {
   for(size_t n = 0; n < blockID.size(); n++)	// This relabelling is the inefficient bit.
    if( blockID[n] == blockID[idx2] )
     blockID[n] = blockID[idx1];
   }
  }
 }
} // namespace Private_PosConstr


template <typename T>
LinearAlgebra::Vector<sdp_gs::Indices> blockStructure(PositivityConstraints<T> const& constr)
 {
 sdp_gs::Indices blockID( constr.dim(), 0 );
 sdp_gs::Index   maxBlock = 0;
 for(size_t n = 0; n < constr.sizeReal(); n++)
  Private_PosConstr::update_blockIDs<T>( blockID, maxBlock, constr.matricesReal()[n] );
 for(size_t n = 0; n < constr.sizeImag(); n++)
  Private_PosConstr::update_blockIDs<T>( blockID, maxBlock, constr.matricesImag()[n] );
 
 LinearAlgebra::Vector<sdp_gs::Indices> ret;
 typedef std::map<sdp_gs::Index,sdp_gs::Index> MapIdxIdx;
 MapIdxIdx blockNumber;
 for(sdp_gs::Index n = 0; n < blockID.size(); n++)
  if( blockID[n] != 0 )		// Ignore completely unused operator basis elements.
   {
   MapIdxIdx::const_iterator i = blockNumber.find( blockID[n] );
   if( i != blockNumber.end() )  ret[i->second] << n;
   else
    {
    blockNumber[ blockID[n] ] = ret.size();
    ret << sdp_gs::Indices( 1, n );
    }
   }
 return ret;
 }


//--------------------------------------------------------------------[ separate_blocks ]-
template <typename T>
LinearAlgebra::Vector<PositivityConstraints<T> >
separate_blocks(PositivityConstraints<T> const& constr, LinearAlgebra::Vector<sdp_gs::Indices> const& blockIdcs)
 {
 using namespace LinearAlgebra;
 Vector<PositivityConstraints<T> > ret( blockIdcs.size(),
                                        PositivityConstraints<T>( constr.sizeReal(), constr.sizeImag() ) );
 for(size_t b = 0; b < blockIdcs.size(); b++)
  {
  for(size_t n = 0; n < constr.sizeReal(); n++)
   ret[b].matricesReal()[n] = get_sparse_block( constr.matricesReal()[n], blockIdcs[b], blockIdcs[b] );
  for(size_t n = 0; n < constr.sizeImag(); n++)
   ret[b].matricesImag()[n] = get_sparse_block( constr.matricesImag()[n], blockIdcs[b], blockIdcs[b] );
  }
 return ret;
 }


//----------------------------------------------------------------[ GreenBasisStructure ]-
namespace Private_PosConstr
{
template <typename T>
LinearAlgebra::Vector<sdp_gs::Indices>
GreenBasisStructure(LinearAlgebra::Vector<typename PositivityConstraints<T>::data_type const*> const& constr,
                    LinearAlgebra::Vector<sdp_gs::Indices>& constraintIdcs)
 {
 constraintIdcs.resize( 0 );
 if( !constr.size() )  return LinearAlgebra::Vector<sdp_gs::Indices>();
 
 LinearAlgebra::Vector<sdp_gs::Indices> ret;
 LabelList<sdp_gs::Indices> chr;
 
 for(size_t n = 0; n < constr[0]->size(); n++)
  {
  sdp_gs::Indices cIdcs;
  for(size_t b = 0; b < constr.size(); b++)
   if( nnz( (*constr[b])[n] ) )  cIdcs << b;
  size_t ng = chr.insert( cIdcs );
  if( ng >= ret.size() )  ret << sdp_gs::Indices();
  ret[ng] << n;
  }
 constraintIdcs = chr.labels(); 
 CHECK_EQUAL( constraintIdcs.size(), ret.size() );
 return ret;
 }
} // namespace Private_PosConstr


template <typename T>
LinearAlgebra::Vector<sdp_gs::Indices>
GreenBasisStructure_Re(LinearAlgebra::Vector< PositivityConstraints<T> > const& constr,
                       LinearAlgebra::Vector<sdp_gs::Indices>& constraintIdcs)
 {
 LinearAlgebra::Vector<typename PositivityConstraints<T>::data_type const*> c( constr.size() );
 for(size_t b = 0; b < constr.size(); b++)
  {
  CHECK_EQUAL( constr[b].sizeReal(), constr[0].sizeReal() );
  c[b] = &constr[b].matricesReal();
  }
 return Private_PosConstr::GreenBasisStructure<T>( c, constraintIdcs );
 }

template <typename T>
LinearAlgebra::Vector<sdp_gs::Indices>
GreenBasisStructure_Im(LinearAlgebra::Vector< PositivityConstraints<T> > const& constr,
                       LinearAlgebra::Vector<sdp_gs::Indices>& constraintIdcs)
 {
 LinearAlgebra::Vector<typename PositivityConstraints<T>::data_type const*> c( constr.size() );
 for(size_t b = 0; b < constr.size(); b++)
  {
  CHECK_EQUAL( constr[b].sizeImag(), constr[0].sizeImag() );
  c[b] = &constr[b].matricesImag();
  }
 return Private_PosConstr::GreenBasisStructure<T>( c, constraintIdcs );
 }


//##################################################################################################
