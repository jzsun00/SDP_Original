//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 04.07.2011 Thomas Barthel
//% SDP groundstate kit
//% Subsystems {\Omega_k}_k.
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// key words: TODO, QUESTION, TEST, CHECK, REMARK, CAUTION
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#ifndef SDP_GS_KIT_SUBSYSTEMS00_H
#define SDP_GS_KIT_SUBSYSTEMS00_H

#include "basis/QuantumNumber00.h"



//__________________________________________________________________________________________________
//##################################################################################### Subsystems #
// Please note that the subsystem index k < NoSubsystems() is shifted by -1 in comparison to the corresponding
//  publication arXiv:1106.4966.
class Subsystems
 {
 public:
 typedef LinearAlgebra::Vector<std::string> strings_type;
 typedef LinearAlgebra::Vector<Subsystem>   container_type;
 typedef LabelList<std::string>             StringList;
 
 Subsystems()  {}
 Subsystems(strings_type const& i_qnNames, sdp_gs::Indices const& i_labelIdcs, QuantumNumbers const& i_basis);
 Subsystems(strings_type const& i_qnNames, sdp_gs::Indices const& i_labelIdcs, QuantumNumbers const& i_basis,
            container_type const& i_Omega);
 
 Subsystems& operator=(Subsystems const& v)  { Basis = v.Basis;  Omega = v.Omega;
                                               QnNames = v.QnNames;  LabelIdcs = v.LabelIdcs;  return *this; }
 
 size_t size()         const  { return Basis.size(); }
 size_t NoSubsystems() const  { return Omega.size(); }
 OneParticleBasis const& basis()              const  { return Basis; }
 container_type   const& operator()()         const  { return Omega; }
 Subsystem        const& operator[](size_t k) const  { return Omega[k]; }
 StringList       const& qnNames()            const  { return QnNames; }
 sdp_gs::Indices  const& labelIdcs()          const  { return LabelIdcs; }
 QuantumNumbers          modeLabels()         const  { return qn_projected( Basis.labels(), LabelIdcs ); }
 
 bool is_complete() const
  { return( Omega == container_type( NoSubsystems(), Subsystem( LinearAlgebra::range(0,size()) ) ) ); }
 
 void append_subsystem(Subsystem      const& OmegaK);
 void append_subsystem(QuantumNumbers const& OmegaK_qns);
 void append_system()  { append_subsystem( Subsystem( LinearAlgebra::range(0,size()) ) ); }
 void set_labelIdcs(sdp_gs::Indices const& i_labelIdcs)  { LabelIdcs = i_labelIdcs; }
 
 PStream::opstream& write(PStream::opstream &out) const;
 PStream::ipstream& read (PStream::ipstream &in);
 
 protected:
 static void check_sets(size_t NoModes, container_type const& i_Omega);
 
 StringList      QnNames;
 sdp_gs::Indices LabelIdcs;
 
 OneParticleBasis Basis;
 container_type   Omega;
 };


//-----------------------------------------------------------------------[ stream stuff ]-
inline PStream::opstream& operator<<(PStream::opstream& out, Subsystems const& v)  { return v.write(out); }
inline PStream::ipstream& operator>>(PStream::ipstream& in , Subsystems      & v)  { return v.read (in);  }
std::ostream& operator<<(std::ostream& out, Subsystems const& v);



//##################################################################################################

#endif // SDP_GS_KIT_SUBSYSTEMS00_H
