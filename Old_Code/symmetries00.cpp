#include "symmetries00.h"

#include "linalg/TensorL02.h"


//__________________________________________________________________________________________________
//######################################################################################## Lattice #

//-----------------------------------------------------------------------[ stream stuff ]-
template<>
const unsigned sdp_gs::ClassID<Lattice>::ID = 7553711;

std::ostream& operator<<(std::ostream& out, Lattice const& v)
 {
 out << "(basis=" << v.basis() << ", L=" << v.L() << ")";
 return out;
 }

PStream::opstream& Lattice::write(PStream::opstream &out) const
 {
 out << (unsigned) sdp_gs::ClassID<Lattice>::ID;	// ID
 out << (int)0; 					// version
 out << A << LV;
 return out;
 }

PStream::ipstream& Lattice::read(PStream::ipstream &in)
 {
 unsigned ID;
 in >> ID;
 CHECK_EQUAL( ID, sdp_gs::ClassID<Lattice>::ID );
 int ver;
 in >> ver;
 DEBUG_CHECK_EQUAL(ver,0);
 in >> A >> LV;
 return in;
 }


//-----------------------------------------------------------------------------[ modulo ]-
sdp_gs::PointF Lattice::modulo(sdp_gs::PointF const& r) const
 {
 using namespace sdp_gs;
 PointF ret(r);
 DEBUG_CHECK_EQUAL( r.size(), dim() );
 PointF::const_iterator LI = LV.begin();
 for(PointF::iterator rI = ret.begin(); rI != ret.end(); ++rI, ++LI)
  if( *LI )
   *rI = *rI - (*LI)*(int)floor( *rI / *LI );
 return ret;
 }


//-----------------------------------------------------------------------[ sublattice_* ]-
LinearAlgebra::Vector<sdp_gs::PointF>
sublattice_cube(Lattice const& latt, sdp_gs::PointI const& N, sdp_gs::PointF const& x0)
 {
 using namespace LinearAlgebra;
 Vector<sdp_gs::PointF> modes = enumerate( TensorIndex( N ) );
 for(size_t n = 0; n < modes.size(); n++)
  modes[n] = latt.modulo( latt.basis()*modes[n] - x0 );
 return modes;
 }



//__________________________________________________________________________________________________
//######################################################################## common group generators #
namespace sdp_gs
{

Float Cos(Float angle)
 {
 switch(int(4*angle))
  {
  case 0: return 1;
  case 1: return 0;
  case 2: return -1;
  case 3: return 0;
  case 4: return 1;
  }
 Float ret = cos(angle*2*M_PI);
 return ret;
 }

Float Sin(Float angle)
 {
 switch(int(4*angle))
  {
  case 0: return 0;
  case 1: return 1;
  case 2: return 0;
  case 3: return -1;
  case 4: return 0;
  }
 Float ret = sin(angle*2*M_PI);
 return ret;
 }


LinearAlgebra::Matrix<Float> identity(size_t D)
 {
 return LinearAlgebra::ScalarMatrix<int>(D,D,1);
 }

LinearAlgebra::Matrix<Float> rotation_2D(Float angle)
 {
 LinearAlgebra::Matrix<Float> ret(2,2);
 ret(0,0) = ret(1,1) = Cos(angle);
 ret(0,1) = -Sin(angle);
 ret(1,0) =  Sin(angle);
 return ret;
 }

LinearAlgebra::Matrix<Float> rotation_3D(Float angle, PointF const& axis)
 {
 using namespace LinearAlgebra;
 PointF axisN = 1./norm_frob(axis) * axis;
 Matrix<Float> Axis(3,1);
 Axis(all,0) = axisN;
 Matrix<Float> Ucross(3,3,0);
 Ucross(1,0) =  axisN[2];
 Ucross(2,0) = -axisN[1];
 Ucross(2,1) =  axisN[0];
 Ucross -= herm(Ucross);
 return Cos(angle)*identity(3) + Sin(angle)*Ucross + (1-Cos(angle))*(Axis*herm(Axis));
 }

LinearAlgebra::Matrix<Float> reflection(PointF const& normal)
 {
 using namespace LinearAlgebra;
 size_t D = normal.size();
 PointF normalN = 1./norm_frob(normal) * normal;
 Matrix<Float> Normal(D,1);
 Normal(all,0) = normalN;
 return identity(D) - 2*(Normal*herm(Normal));
 }

} // namespace sdp_gs



//__________________________________________________________________________________________________
//######################################################################### common symmetry groups #

//--------------------------------------------------------------[ symmGroup_point_cubic ]-
LinearAlgebra::Vector<sdp_gs::MatrixF> symmGenerators_point_cubic(size_t D)
 {
 using namespace LinearAlgebra;
 CHECK_COMPARE( D, <=, 3 );
 
 Vector<sdp_gs::MatrixF> generators;
 sdp_gs::MatrixF Id = sdp_gs::identity(D);
 
 if( D == 2 )  generators << sdp_gs::rotation_2D( 1./4 );
 if( D == 3 ) 
  for(size_t d = 0; d < D; d++)
   generators << sdp_gs::rotation_3D( 1./4, Id(all,d) );
 
 for(size_t d = 0; d < D; d++)
  generators << sdp_gs::reflection( Id(all,d) );

 return generators;
 }


//-----------------------------------------------------------------[ symmGroup_spinFlip ]-
LinearAlgebra::Vector<sdp_gs::MatrixF> symmGenerators_spinFlip()
 {
 return LinearAlgebra::Vector<sdp_gs::MatrixF>( 1, sdp_gs::reflection( sdp_gs::PointF(1,1) ) );
 }



//__________________________________________________________________________________________________
//##################################################################################### Symmetries #

//------------------------------------------------------------------------[ constructor ]-
Symmetries::Symmetries(bool i_realF, lattice_type const& i_latt, groups_type const& i_groups,
                       sdp_gs::Indices const& i_conservedQnIdcs)
  : realF(i_realF), Latt(i_latt), Groups(i_groups), ConservedQnIdcs(i_conservedQnIdcs)
 {
 #ifndef SDP_GS_COMPLEX
 PRECONDITION( realF == 1 );
 #endif
 CHECK_EQUAL( i_latt.idcs.size(), i_latt.symm.dim() );
 for(size_t n = 0; n < i_groups.size(); n++)
  { CHECK_EQUAL( i_groups[n].idcs.size(), i_groups[n].symm.dim() ); }
 }


//-----------------------------------------------------------------------[ stream stuff ]-
template<>
const unsigned sdp_gs::ClassID<Symmetries>::ID = 7553713;

std::ostream& operator<<(std::ostream& out, Symmetries const& v)
 {
 out << "(is_real="          << v.is_real() << ",\n"
     << " lattice="         << v.lattice() << ",\n"
     << " finite groups="   << v.groups()  << ",\n"
     << " conservedQnIdcs=" << v.conservedQnIdcs()
     << ")";
 return out;
 }

PStream::opstream& Symmetries::write(PStream::opstream &out) const
 {
 out << (unsigned) sdp_gs::ClassID<Symmetries>::ID;	// ID
 out << (int)1; 					// version
 out << realF << Latt << Groups << ConservedQnIdcs;
 return out;
 }

PStream::ipstream& Symmetries::read(PStream::ipstream &in)
 {
 unsigned ID;
 in >> ID;
 CHECK_EQUAL( ID, sdp_gs::ClassID<Symmetries>::ID );
 int ver;
 in >> ver;
 DEBUG_CHECK_EQUAL(ver,1);
 in >> realF >> Latt >> Groups >> ConservedQnIdcs;
 return in;
 }


//-------------------------------------------------------------------------[ operator== ]-
bool Symmetries::operator==(Symmetries const& v) const
 {
 if( Groups.size() != v.Groups.size() )  return false;
 if( !( realF == v.realF && Latt == v.Latt && ConservedQnIdcs == v.ConservedQnIdcs ) )  return false;
 for(size_t n = 0; n < Groups.size(); n++)
  if( !(Groups[n] == v.Groups[n]) )  return false;
 return true;
 }


//--------------------------------------------------------------------[ NoGroupElements ]-
sdp_gs::Indices Symmetries::NoGroupElements() const
 {
 sdp_gs::Indices ret( NoGroups() );
 for(size_t ng = 0; ng < NoGroups(); ng++)  ret[ng] = group(ng).symm.order();
 return ret;
 }


//##################################################################################################
