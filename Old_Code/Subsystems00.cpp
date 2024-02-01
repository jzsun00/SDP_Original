#include "Subsystems00.h"


//__________________________________________________________________________________________________
//##################################################################################### Subsystems #

//-----------------------------------------------------------------------[ stream stuff ]-
template<>
const unsigned sdp_gs::ClassID<Subsystems>::ID = 7553701;

std::ostream& operator<<(std::ostream& out, Subsystems const& v)
 {
 out << "(qnNames="    << v.qnNames() << ",\n"
     << " lableNames=" << v.qnNames()[v.labelIdcs()] << ",\n"
     << " basis="      << v.basis()   << ",\n";
 for(size_t k = 0; k < v.NoSubsystems(); k++)
  out << " Omega_" << k << "=" << v.modeLabels()[ v[k] ] << (k+1<v.NoSubsystems()? ",\n":"");
 out << ")";
 return out;
 }

PStream::opstream& Subsystems::write(PStream::opstream &out) const
 {
 out << (unsigned) sdp_gs::ClassID<Subsystems>::ID;	// ID
 out << (int)1; 					// version
 out << QnNames << LabelIdcs << Basis << Omega;
 return out;
 }

PStream::ipstream& Subsystems::read(PStream::ipstream &in)
 {
 unsigned ID;
 in >> ID;
 CHECK_EQUAL( ID, sdp_gs::ClassID<Subsystems>::ID );
 int ver;
 in >> ver;
 DEBUG_CHECK_EQUAL(ver,1);
 in >> QnNames >> LabelIdcs >> Basis >> Omega;
 return in;
 }


//------------------------------------------------------------------------[ constructor ]-
Subsystems::
Subsystems(strings_type const& i_qnNames, sdp_gs::Indices const& i_labelIdcs, QuantumNumbers const& i_basis)
 : QnNames(i_qnNames), LabelIdcs(i_labelIdcs), Basis(i_basis)
 {
 if( !LabelIdcs.size() )  LabelIdcs = LinearAlgebra::range( 0, i_qnNames.size() );
 if( i_basis.size() )  { CHECK_EQUAL( i_basis[0].size(), i_qnNames.size() ); }
 PRECONDITION( LinearAlgebra::is_set( i_basis ) );
 }

Subsystems::
Subsystems(strings_type const& i_qnNames, sdp_gs::Indices const& i_labelIdcs, QuantumNumbers const& i_basis,
                       container_type const& i_Omega)
 : QnNames(i_qnNames), LabelIdcs(i_labelIdcs), Basis(i_basis), Omega(i_Omega)
 {
 if( !LabelIdcs.size() )  LabelIdcs = LinearAlgebra::range( 0, i_qnNames.size() );
 if( i_basis.size() )  { CHECK_EQUAL( i_basis[0].size(), i_qnNames.size() ); }
 PRECONDITION( LinearAlgebra::is_set( i_basis ) );
 check_sets( i_basis.size(), i_Omega );
 }


//-------------------------------------------------------------------------[ check_sets ]-
void Subsystems::check_sets(size_t NoModes, container_type const& i_Omega)
 {
 using namespace LinearAlgebra;
 if( i_Omega.size() )
  {
  PRECONDITION( is_set( i_Omega[0] ) );
  PRECONDITION( is_subset( i_Omega[0], sdp_gs::Indices(range(0,NoModes)) ) );
  for(size_t k = 1; k < i_Omega.size(); k++)
   {
   PRECONDITION( is_set( i_Omega[k] ) );
   PRECONDITION( is_subset( i_Omega[k], i_Omega[k-1] ) );
   }
  }
 }


//-------------------------------------------------------------------[ append_subsystem ]-
void Subsystems::append_subsystem(Subsystem const& OmegaK)
 {
 using namespace LinearAlgebra;
 PRECONDITION( is_set( OmegaK ) );
 /*if( NoSubsystems() )
  { PRECONDITION( is_subset( OmegaK, Omega[NoSubsystems()-1] ) ); }
 else
  { PRECONDITION( is_subset( OmegaK, sdp_gs::Indices(range(0,size())) ) ); }*/
 PRECONDITION( is_subset( sortC(OmegaK), sdp_gs::Indices(range(0,size())) ) );
 Omega << OmegaK;
 }

void Subsystems::append_subsystem(QuantumNumbers const& OmegaK_qns)
 {
 using namespace LinearAlgebra;
 PRECONDITION( is_set( OmegaK_qns ) );
 PRECONDITION( Basis.contains( OmegaK_qns ) );
 append_subsystem( Basis( OmegaK_qns ) );
 }


//##################################################################################################
