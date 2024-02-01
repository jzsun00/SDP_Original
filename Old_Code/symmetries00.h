//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//% 07.07.2011 Thomas Barthel
//% SDP groundstate kit
//% Symmetries of the system (translation invariance, point sysmmetries)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// key words: TODO, QUESTION, TEST, CHECK, REMARK, CAUTION
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#ifndef SDP_GS_KIT_SYMMETRIES00_H
#define SDP_GS_KIT_SYMMETRIES00_H

#include "configuration.h"

#include "linearalgebra/scalarmatrix.h"
#include "linearalgebra/pstreamio.h"



//##################################################################################################

//-----------------------------------------------------------------------[ distance_min ]-
template <typename T, typename L>
double distance_min(T const& x, L const& list,
  typename boost::enable_if< LinearAlgebra::is_vector<L> >::type* dummy1=0,
  typename boost::enable_if< boost::is_same<T, typename LinearAlgebra::interface<L>::value_type> >::type* dummy2=0);

template <typename T, typename L>
double distance_min(T const& x, L const& list, size_t& closest_element,
  typename boost::enable_if< LinearAlgebra::is_vector<L> >::type* dummy1=0,
  typename boost::enable_if< boost::is_same<T, typename LinearAlgebra::interface<L>::value_type> >::type* dummy2=0);

template <typename L1, typename L2>
double distance_min(L1 const& list1, L2 const& list2,
  typename boost::enable_if< LinearAlgebra::is_vector<L1> >::type* dummy1=0,
  typename boost::enable_if< LinearAlgebra::is_vector<L2> >::type* dummy2=0,
  typename boost::enable_if< boost::is_same<typename LinearAlgebra::interface<L1>::value_type,
                                            typename LinearAlgebra::interface<L2>::value_type> >::type* dummy3=0);

template <typename L>
double distance_min(L const& list,
  typename boost::enable_if< LinearAlgebra::is_vector<L> >::type* dummy1=0);



//__________________________________________________________________________________________________
//######################################################################################## Lattice #
class Lattice
 {
 public:
 typedef sdp_gs::MatrixF basis_type;
 
 Lattice()  {}
 explicit
 Lattice(sdp_gs::PointF const& i_LV)
  : A(LinearAlgebra::ScalarMatrix<int>(i_LV.size(),i_LV.size(),1)), LV(i_LV)  {}
 Lattice(basis_type const& i_A, sdp_gs::PointF const& i_LV) : A(i_A), LV(i_LV)
  { CHECK_EQUAL( i_A.size1(), i_LV.size() ); }
 bool operator==(Lattice const& v) const  { return( A == v.A && LV == v.LV ); }  // "!(==)" is not (yet) "!="
 
 size_t                dim()   const  { return LV.size(); }
 basis_type const&     basis() const  { return A; }
 sdp_gs::PointF const& L()     const  { return LV; }
 
 sdp_gs::PointF modulo(sdp_gs::PointF const& r) const;
 
 PStream::opstream& write(PStream::opstream &out) const;
 PStream::ipstream& read (PStream::ipstream &in);
 
 protected:
 basis_type	A;	// basis vectors of the lattice (only for the directions for which we have periodicity)
 sdp_gs::PointF	LV;	// entries = 0 indicate 'no PBC' for that (Euclidian) direction
 };


//-----------------------------------------------------------------------[ stream stuff ]-
inline PStream::opstream& operator<<(PStream::opstream& out, Lattice const& v)  { return v.write(out); }
inline PStream::ipstream& operator>>(PStream::ipstream& in , Lattice      & v)  { return v.read (in);  }
std::ostream& operator<<(std::ostream& out, Lattice const& v);


//-----------------------------------------------------------------------[ sublattice_* ]-
LinearAlgebra::Vector<sdp_gs::PointF>
sublattice_cube(Lattice const& latt, sdp_gs::PointI const& N, sdp_gs::PointF const& x0);

inline LinearAlgebra::Vector<sdp_gs::PointF>
sublattice_cube(Lattice const& latt, sdp_gs::PointI const& N)
 { return sublattice_cube( latt, N, sdp_gs::PointF( latt.dim(), 0 ) ); }



//__________________________________________________________________________________________________
//######################################################################## common group generators #
namespace sdp_gs
{
sdp_gs::MatrixF identity   (size_t D);
sdp_gs::MatrixF rotation_2D(sdp_gs::Float angle);				// angle in units of 2pi
sdp_gs::MatrixF rotation_3D(sdp_gs::Float angle, sdp_gs::PointF const& axis);	// angle in units of 2pi
sdp_gs::MatrixF reflection (sdp_gs::PointF const& normal);
} // namespace sdp_gs



//__________________________________________________________________________________________________
//####################################################################################### GroupRep #
// A class holding the linear representation of a symmetry group, as applicable to (some subset)
//  of the mode quantum numbers.
template <typename element_type>
class GroupRep
 {
 public:
 typedef LinearAlgebra::Vector<element_type> container_type;
 
 GroupRep()  {}
 template <typename T>
 GroupRep(T const& i_elements, typename boost::enable_if<LinearAlgebra::is_vector<T> >::type* dummy = 0);
 bool operator==(GroupRep<element_type> const& v) const  // "!(==)" is not (yet) "!="
  { return( Elements == v.Elements ); }
 
 size_t order() const  { return Elements.size(); }
 size_t dim()   const  { return order()? Elements[0].size1() : 0; }
 
 element_type   const& operator[](size_t n) const  { return Elements[n]; }
 container_type const& elements()           const  { return Elements; }
 
 template <typename T>
 inline void generate(T const& generators, element_type const& Id,
                      typename boost::enable_if<LinearAlgebra::is_vector<T> >::type* dummy = 0);
 
 PStream::opstream& write(PStream::opstream &out) const  { out << Elements;  return out; }
 PStream::ipstream& read (PStream::ipstream &in)         { in  >> Elements;  return in; }
 
 protected:
 container_type	Elements;
 };


//-----------------------------------------------------------------------[ stream stuff ]-
template <typename E>
inline PStream::opstream& operator<<(PStream::opstream& out, GroupRep<E> const& v)  { return v.write(out); }
template <typename E>
inline PStream::ipstream& operator>>(PStream::ipstream& in , GroupRep<E>      & v)  { return v.read (in);  }
template <typename E>
inline std::ostream& operator<<(std::ostream& out, GroupRep<E> const& v)
 { out << small2zero( v.elements() );  return out; }



//__________________________________________________________________________________________________
//######################################################################### common symmetry groups #
inline GroupRep<sdp_gs::MatrixF> generate_group(LinearAlgebra::Vector<sdp_gs::MatrixF> const& generators,
                                                sdp_gs::MatrixF const& Id);
inline GroupRep<sdp_gs::MatrixF> generate_group(LinearAlgebra::Vector<sdp_gs::MatrixF> const& generators, size_t Dim)
 { return generate_group( generators, sdp_gs::identity(Dim) ); }

LinearAlgebra::Vector<sdp_gs::MatrixF> symmGenerators_point_cubic(size_t D);
LinearAlgebra::Vector<sdp_gs::MatrixF> symmGenerators_spinFlip();

inline GroupRep<sdp_gs::MatrixF> symmGroup_point_cubic(size_t D)
 { return generate_group( symmGenerators_point_cubic(D), D ); }
inline GroupRep<sdp_gs::MatrixF> symmGroup_spinFlip()
 { return generate_group( symmGenerators_spinFlip(), 1 ); }



//__________________________________________________________________________________________________
//################################################################################## IndexSymmetry #
// A struct, that tells us to which quantum numbers, the corresponding symmetry applies.
template <typename Symm>
struct IndexSymmetry
 {
 IndexSymmetry()  {}
 IndexSymmetry(sdp_gs::Indices const& i_idcs, Symm const& i_symm) : idcs(i_idcs), symm(i_symm)
  { PRECONDITION( LinearAlgebra::is_set( i_idcs ) ); }
 bool operator==(IndexSymmetry<Symm> const& v) const  // "!(==)" is not (yet) "!="
  { return( idcs == v.idcs && symm == v.symm ); }
 
 sdp_gs::Indices idcs;
 Symm            symm;
 };


//-----------------------------------------------------------------------[ stream stuff ]-
template <typename S>
inline PStream::opstream& operator<<(PStream::opstream& out, IndexSymmetry<S> const& v)
 { out << v.idcs << v.symm;  return out; }
template <typename S>
inline PStream::ipstream& operator>>(PStream::ipstream& in , IndexSymmetry<S>      & v)
 { in  >> v.idcs >> v.symm;  return in;  }
template <typename S>
inline std::ostream&      operator<<(std::ostream& out,      IndexSymmetry<S> const& v)
 { out << "(indices=" << v.idcs << ", symmetry=" << v.symm << ")";  return out; }



//__________________________________________________________________________________________________
//##################################################################################### Symmetries #
// This class holds all symmetries corresponding to a quantum lattice system in second quantization:
//  translation symmetries, point symmetries, and conserved quantum numbers.
// The (quantum number) indices refer to those used in the Subsystems class.
// For the moment, we assume that elements of different groups commute!
class Symmetries
 {
 public:
 typedef IndexSymmetry<Lattice>            lattice_type;
 typedef sdp_gs::MatrixF                   operator_type;
 typedef GroupRep<operator_type>           groupRep_type;
 typedef IndexSymmetry<groupRep_type>      group_type;
 typedef LinearAlgebra::Vector<group_type> groups_type;
 
 Symmetries()  {}
 Symmetries(bool i_realF, lattice_type const& i_latt, groups_type const& i_groups,
            sdp_gs::Indices const& i_conservedQnIdcs);
 bool operator==(Symmetries const& v) const;  // "!(==)" is not (yet) "!="
 
 lattice_type    const& lattice()         const  { return Latt; }
 size_t                 NoGroups()        const  { return Groups.size(); }
 sdp_gs::Indices        NoGroupElements() const;
 size_t                 NoGroupElementsTot() const  { return LinearAlgebra::prod( NoGroupElements() ); }
 group_type      const& group(size_t n)   const  { return Groups[n]; }
 groups_type     const& groups()          const  { return Groups; }
 sdp_gs::Indices const& conservedQnIdcs() const  { return ConservedQnIdcs; }
 
 bool is_real() const          { return realF; }
 void set_realF(bool i_realF)  { realF = i_realF; }
 
 void append_group(group_type const& group)  { Groups << group; }
 
 PStream::opstream& write(PStream::opstream &out) const;
 PStream::ipstream& read (PStream::ipstream &in);
 
 protected:
 bool            realF;		// Hamiltonian is real in the chosen basis => Im G_s = 0 for all s.
 lattice_type    Latt;
 groups_type     Groups;
 sdp_gs::Indices ConservedQnIdcs;
 };


//-----------------------------------------------------------------------[ stream stuff ]-
inline PStream::opstream& operator<<(PStream::opstream& out, Symmetries const& v)  { return v.write(out); }
inline PStream::ipstream& operator>>(PStream::ipstream& in , Symmetries      & v)  { return v.read (in);  }
std::ostream& operator<<(std::ostream& out, Symmetries const& v);



//##################################################################################################


#include "symmetries00.cc"

#endif // SDP_GS_KIT_SYMMETRIES00_H
