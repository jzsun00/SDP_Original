#include "algebra/operator00.h"
#include "latticeAlgebra/operatorTransform00.h"
#include "io/symmetries00_io.h"
#include "io/Subsystems00_io.h"

#include "fileioL00.h"
#include "parseMathL03.h"
#include "labelledValueL02.h"
#include "linalg/TensorL02.h"


using namespace std;
using namespace LinearAlgebra;



//##################################################################################################

//----------------------------------------------------------------------------[ ModeSet ]-
struct ModeSet
 {
 ModeSet() : copyF(0), k(0), cond(1), orderFct(0), relativeVolume(-1)  {}
 
 void set_copyFlag()  { copyF = 1+2+4+8+16; }
 
 void append_N(ParseMath::Formula const& f)
  { if( copyF&1 ) { copyF^=1; N.resize(0); }  N << f; }
 void append_trafo_n(ParseMath::Formula const& f)
  { if( copyF&2 ) { copyF^=2; trafo_n.resize(0); }  trafo_n << f; }
 void append_trafo_coord(ParseMath::Formula const& f)
  { if( copyF&4 ) { copyF^=4; trafo_coord.resize(0); }  trafo_coord << f; }
 void append_qnName(string const& f)
  { if( copyF&8 ) { copyF^=8; qnName.resize(0); }  qnName << f; }
 void append_qnValue(ParseMath::Formula const& f)
  { if( copyF&16) { copyF^=16; qnValue.resize(0); }  qnValue << f; }
 
 long   copyF;
 size_t k;
 Vector<ParseMath::Formula> N;
 Vector<ParseMath::Formula> trafo_n;
 Vector<ParseMath::Formula> trafo_coord;
 Vector<string>             qnName;
 Vector<ParseMath::Formula> qnValue;
 ParseMath::Formula         cond;
 ParseMath::Formula         orderFct;
 ParseMath::Formula         relativeVolume;
 };

std::ostream& operator<<(std::ostream& out, ModeSet const& v)
 {
 out << "(k="  << v.k
     << ", N=" << v.N
     << ", trafo_n="     << v.trafo_n
     << ", trafo_coord=" << v.trafo_coord
     << ", qnName="  << v.qnName
     << ", qnValue=" << v.qnValue
     << ", cond="    << v.cond
     << ", orderFct="<< v.orderFct
     << ", relativeVolume=" << v.relativeVolume
     << ")";
 return out;
 }


//---------------------------------------------------------------[ varMap_insert_vector ]-
template <typename T>
void varMap_insert_vector(ParseMath::VariableMap& varMap, T const& v, string const& label)
 {
 for(size_t i = 0; i < v.size(); i++)
  varMap[ label + boost::lexical_cast<string>(i) ] = v[i];
 }

template <typename T>
void varMap_insert_vector(ParseMath::VariableMap& varMap, T const& v, Vector<string> const& label)
 {
 CHECK_EQUAL( v.size(), label.size() );
 for(size_t i = 0; i < v.size(); i++)
  varMap[ label[i] ] = v[i];
 }


//--------------------------------------------------------------------------[ sort_symm ]-
// Sorting modes according to some order function and trying to put modes of the same symmetry orbit
//  and order function value close to each other.
typedef pair<QuantumNumber,double> OrderedMode;
class CompareOrderedMode
 {
 public:
 bool operator()(OrderedMode const& a, OrderedMode const& b)
  { return( a.second < b.second || ( a.second == b.second && a.first < b.first ) ); }
 };

QuantumNumbers sort_symm(Symmetries const& symm, QuantumNumbers const& qn, Vector<double> const& order)
 {
 CHECK_EQUAL( qn.size(), order.size() );
 typedef std::set<OrderedMode,CompareOrderedMode> OrderedModeSet;
 OrderedModeSet modes;
 for(size_t n = 0; n < qn.size(); n++)  modes.insert( OrderedMode( qn[n], order[n] ) );
 
 QuantumNumbers ret( qn.size() );
 size_t m = 0;
 while( modes.size() )
  {
  OrderedMode mode = *modes.begin();
  modes.erase( modes.begin() );
  ret[m++] = mode.first;
  QuantumNumber rep = representative( mode.first, symm );
  for(OrderedModeSet::iterator i = modes.begin(); i != modes.end() && i->second == mode.second; ++i)
   if( representative( i->first, symm ) == rep )
    {
    ret[m++] = i->first;
    modes.erase( i );
    }
  }
 CHECK_EQUAL( m, qn.size() );
 return ret;
 }



//__________________________________________________________________________________________________
//########################################################################################### MAIN #
int main(int argc, char* argv[])
 {
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 
 string program = "create-subsystems";
 
 cout<<" >>> Create Subsystems object <<< 11.08.2011\n";
 cout<<" usage:  "<< program << endl;
 cout<<"		 -in_symm <*.Symm.bin file, providing Lattice and qnNames>  (-show_history)  (-show_params)\n";
 cout<<"		(-extra_param <param label> <formula> ...)\n";
 cout<<"		(-label_subsys <label for the Subsystems object>)  (-dont_append_paramExt_to_label)\n";
 cout<<"		(-volume <k=monomOrder-1> <formula for total number of sites =-1> ...)\n";
 cout<<"		 -add_modes  (-copy_modes <index to the mode set =-1>)\n";
 cout<<"		(-m_k <monomOrder-1 =0>)  (-m_cubeSize <formulas for cube sizes N0=11 ...)>)\n";
 cout<<"		(-m_cubeTrafo <formula for n0' ...>)  (-m_coordTrafo <formula for coord0'> ...)  (-m_cond <formula for condition>)\n";
 cout<<"		(-m_qnValue <qnName> <formula> ...)\n";
 cout<<"		(-m_orderFct <formula>)  (-m_relativeVolume <formula for ratio =-1>)\n";
 cout<<"		(-Nadd_fullSystem <number of times to add the full system in the end =0>)\n";
 cout<<"		(-qnLabelNames <qnName to be used for labelling> ...)  (-dont_show_subsystems)\n";
 cout<<"		(-out_subsys <*.Subsys.bin file =std string>)\n";
 cout<<" You can add sets of modes to some (monom order = k+1) subsystem by calling -add_modes\n";
 cout<<"  (several times). Each time you call -add_modes, you can/should add corresponding -m_* specifications.\n";
 cout<<"  'formula' can be an arbitrary function of the lattice site numbers n0,n1,.., all quantum numbers,\n";
 cout<<"   of the parameters from the history of the *.Symm.bin file, and of the extra parameters.\n";
 cout<<"  The further quantum numbers (-m_qnValue), mode condition (-m_cond), and ordering function (-m_orderFct)\n";
 cout<<"   are evaluated after applying the cube and coordinate transformations.\n";
 cout<<"  The ordering function (-m_orderFct) and relative volume (-m_relativeVolume) can be used to adjust the number of sites for a mode set\n";
 cout<<"   relative to the requested total number of sites (-volume) for the given monom order (k+1).\n";
 cout<<" If -qnLabelIdcs, is not specified, the real-space coordinates are used.\n";
 cout<<endl;
 
 string in_symm;
 string out_subsys;
 string label_subsys;
 bool   show_history = 0;
 bool   show_params  = 0;
 bool   append_paramExt_to_label = 1;
 bool   show_subsystems = 1;
 size_t Nadd_fullSystem = 0;
 Vector<string>		qnLabelNames;
 Vector<string>		paramExt_label;
 Vector<ParseMath::Formula> paramExt;
 Vector<ModeSet>	modes;
 map<size_t,ParseMath::Formula> volumeF;
 
 int nm = -1;
 
 
 //=================================================================================
 //-[read in command line parameters]-----------------------------------------------
 if( argc<2) return 1;
 
 int argMode = 0;
 set<string> argModable;
 map<string,int> argID;
 argID["-in_symm"]		= 1;
 argID["-show_history"]		= 10;
 argID["-show_params"]		= 12;
 argID["-extra_param"]		= 14;		argModable.insert("-extra_param");
 argID["-label_subsys"]		= 20;
 argID["-dont_append_paramExt_to_label"] = 22;
 argID["-volume"]		= 40;		argModable.insert("-volume");
 argID["-add_modes"]		= 42;
 argID["-copy_modes"]		= 44;
 argID["-m_k"]			= 48;
 argID["-m_cubeSize"]		= 50;		argModable.insert("-m_cubeSize");
 argID["-m_cubeTrafo"]		= 52;		argModable.insert("-m_cubeTrafo");
 argID["-m_coordTrafo"]		= 54;		argModable.insert("-m_coordTrafo");
 argID["-m_qnValue"]		= 56;		argModable.insert("-m_qnValue");
 argID["-m_cond"]		= 58;
 argID["-m_orderFct"]		= 60;
 argID["-m_relativeVolume"]	= 62;
 argID["-Nadd_fullSystem"]	= 70;
 argID["-qnLabelNames"]		= 74;		argModable.insert("-qnLabelNames");
 argID["-dont_show_subsystems"]	= 76;
 argID["-out_subsys"]		= 80;
 
 for(int argn=1; argn<argc; argn++)
  {
  int aID;  string argStr( argv[argn] );
  if( argStr[0] == '-' && (argStr[1]<'0'||argStr[1]>'9') )
   { aID = argID[argStr];  argMode = 0;  if(argModable.count(argStr)) argMode = aID; }
  else { aID=argMode;  argn--; }
  switch( aID )
   {
   case   1:  in_symm = argv[argn+1];  argn++;  break;
   case  10:  show_history = 1;  break;
   case  12:  show_params  = 1;  break;
   case  14:  paramExt_label << argv[argn+1];  argn++;
              paramExt << ParseMath::parse( argv[argn+1] );  argn++;  break;
   case  20:  label_subsys = argv[argn+1];  argn++;  break;
   case  22:  append_paramExt_to_label = 0;  break;
   case  40:  volumeF[ boost::lexical_cast<size_t>( argv[argn+1] ) ] = ParseMath::parse( argv[argn+2] );
              argn+=2;  break;
   case  42:  nm++;  resize( modes, nm+1, true );  break;
   case  44:  nm++;  resize( modes, nm+1, true );
              modes[nm] = modes[ (boost::lexical_cast<long>(argv[argn+1])+nm) % nm ];  argn++;
              modes[nm].set_copyFlag();  break;
   case  48:  modes[nm].k = boost::lexical_cast<size_t>( argv[argn+1] );  argn++;  break;
   case  50:  modes[nm].append_N(           ParseMath::parse( argv[argn+1] ) );  argn++;  break;
   case  52:  modes[nm].append_trafo_n(     ParseMath::parse( argv[argn+1] ) );  argn++;  break;
   case  54:  modes[nm].append_trafo_coord( ParseMath::parse( argv[argn+1] ) );  argn++;  break;
   case  56:  modes[nm].append_qnName(  argv[argn+1] );  argn++;
              modes[nm].append_qnValue( ParseMath::parse( argv[argn+1] ) );  argn++;  break;
   case  58:  modes[nm].cond     = ParseMath::parse( argv[argn+1] );  argn++;  break;
   case  60:  modes[nm].orderFct = ParseMath::parse( argv[argn+1] );  argn++;  break;
   case  62:  modes[nm].relativeVolume = ParseMath::parse( argv[argn+1] );  argn++;  break;
   case  70:  Nadd_fullSystem = boost::lexical_cast<size_t>( argv[argn+1] );  argn++;  break;
   case  74:  qnLabelNames << argv[argn+1];  argn++;  break;
   case  76:  show_subsystems = 0;  break;
   case  80:  out_subsys = argv[argn+1];  argn++;  break;
   default :  PANIC("Wrong command line arguments!")(argStr);  break;
   }
  }
 
 PRECONDITION( fileExists(in_symm) )( in_symm );
 
 
 cout << "\n========================================================== read Symmetries file ==\n";
 History hist;
 LabelList<std::string> qnNames;
 Symmetries      symm;
 Lattice         latt;
 sdp_gs::Indices coord_idcs;
 {
 read( in_symm, hist, symm, qnNames );
 latt       = symm.lattice().symm;
 coord_idcs = symm.lattice().idcs;
 }
 string label_latt = hist[-1].get_info_s("label_Lattice");
 if( show_history )  { TRACE( hist ); }
 TRACE( qnNames )( qnNames[coord_idcs] );
 
 // get parameters from History
 map<string,double> param;
 map<string,string> paramString;	// not used
 hist.export_all( paramString, param );
 ParseMath::VariableMap varMap( param );
 
 // add extra parameters
 for(size_t n = 0; n < paramExt.size(); n++)
  varMap[ paramExt_label[n] ] = ParseMath::eval( paramExt[n], varMap );
 if( show_params )  { TRACE( varMap );  return -1; }
 
 
 cout << "\n======================================================================= History ==\n";
 History::Action a;
 a.info_s["program"]      = program;
 a.info_s["command_line"] = commandLine( argc, argv );
 a.info_s["cwd"]          = current_workingPath();
 a.info_s["date"]         = dateTimeString();
 a.info_d["algebraF"]     = hist[-1].get_info_d("algebraF");
 hist.push_back(a);
 
 for(size_t n = 0; n < paramExt.size(); n++)
  {
  double val = ParseMath::eval_double( paramExt[n], varMap );
  hist[-1].param_d[ paramExt_label[n] ] = val;
  if( append_paramExt_to_label )
   label_subsys += paramExt_label[n] + boost::lexical_cast<string>( val );
  }
 
 hist[-1].info_s["label_Subsystems"] = label_subsys;
 
 
 cout << "\n============================================================= create Subsystems ==\n";
 // 1) Generate all mode reservoirs (not yet balancing the volumes).
 Vector<Vector<QuantumNumbers> > reservoirs;
 Vector<Vector<ModeSet> >        modeSets;
 
 for(size_t m = 0; m < modes.size(); m++)
  {
  TRACE( modes[m] );
  ModeSet const& M( modes[m] );
  if( M.k >= reservoirs.size() )
   {
   resize( reservoirs, M.k+1, true );
   resize( modeSets,   M.k+1, true );
   }
  OneParticleBasis Qn;
  Vector<double>   orderVal;
  
  PRECONDITION( !M.N.size() || M.N.size() == latt.basis().size2() )( M.N );
  Vector<int> N( latt.basis().size2(), 11 );
  for(size_t i = 0; i < M.N.size(); i++)  N[i] = (int)ParseMath::eval_double( M.N[i], varMap );
  TensorIndex ti(N);
  
  for(TensorIndexState n = ti.begin(); n.is_valid(); ++n)
   {
   ParseMath::VariableMap varMap2( varMap );
   
   // transform lattice site index
   Vector<int> nV = n();
   if( M.trafo_n.size() )
    {
    CHECK_EQUAL( M.trafo_n.size(), nV.size() );
    varMap_insert_vector( varMap2, nV, "n" );
    for(size_t i = 0; i < nV.size(); i++)
     nV[i] = (int)ParseMath::eval_double( M.trafo_n[i], varMap2 );
    }
   varMap_insert_vector( varMap2, nV, "n" );
   
   // transform lattice site coordinate
   Vector<sdp_gs::Float> r = latt.basis()*nV;
   if( M.trafo_coord.size() )
    {
    CHECK_EQUAL( M.trafo_coord.size(), r.size() );
    varMap_insert_vector( varMap2, r, qnNames[coord_idcs] );
    for(size_t i = 0; i < r.size(); i++)
     r[i] = ParseMath::eval_double( M.trafo_coord[i], varMap2 );
    r = latt.modulo( r );
    }
   varMap_insert_vector( varMap2, r, qnNames[coord_idcs] );
   
   // add the remaining quantum numbers to varMap
   for(size_t i = 0; i < qnNames.size(); i++)
    if( distance_min( (sdp_gs::Index)i, coord_idcs ) > zeroThresh )
     varMap2[ qnNames[i] ] = 0;
   CHECK_COMPARE( M.qnName.size() + coord_idcs.size(), <=, qnNames.size() );
   for(size_t i = 0; i < M.qnName.size(); i++)
    varMap2[ M.qnName[i] ] = ParseMath::eval_double( M.qnValue[i], varMap2 );
   
   // add mode, if the given condition is true.
   if( ParseMath::eval_double( M.cond, varMap2 ) != 0 )
    {
    QuantumNumber qn( qnNames.size() );
    for(size_t i = 0; i < qnNames.size(); i++)
     qn[i] = ParseMath::eval_double( qnNames[i], varMap2 );
    size_t idx = Qn.insert( qn );
    if( idx == orderVal.size() )  orderVal << ParseMath::eval_double( M.orderFct, varMap2 );
    }
   }
  
  //reservoirs[M.k] << Qn.labels()[ sortIdx( orderVal ) ];
  reservoirs[M.k] << sort_symm( symm, Qn.labels(), orderVal );
  modeSets[M.k]   << M;
  } // next ModeSet
 
 // 2) Generate all mode sets, balancing the volumes, if required.
 Vector<long> volume( reservoirs.size(), -1 );
 for(size_t k = 0; k < reservoirs.size(); k++)
  if( volumeF.find(k) != volumeF.end() )
   volume[k] = (long)ParseMath::eval_double( volumeF[k], varMap );
 
 Vector<OneParticleBasis> Qn( reservoirs.size() );
 for(size_t k = 0; k < reservoirs.size(); k++)
  {
  sdp_gs::Indices adjustable;
  Vector<double> relativeVolume;
  for(size_t n = 0; n < modeSets[k].size(); n++)
   {
   double relVol = ParseMath::eval_double( modeSets[k][n].relativeVolume, varMap );
   if( relVol >= 0 )
    {
    adjustable << n;
    relativeVolume << relVol;
    }
   else Qn[k].insert( reservoirs[k][n] );
   }
  
  TRACE(k)(adjustable)(volume[k]);
  if( adjustable.size() )
   {
   PRECONDITION( volume[k] != -1 );
   CHECK_COMPARE( volume[k], >=, Qn[k].size() );
   CHECK_COMPARE( fabs( sum(relativeVolume) - 1 ), <, zeroThresh );
   sdp_gs::Indices NoSites( adjustable.size(), 0 );
   while( Qn[k].size() < volume[k] )
    {
    Vector<double> ratio( NoSites );
    if( sum(NoSites) )  ratio *= 1./sum(NoSites);
    ratio -= relativeVolume;
    size_t a = min_element( ratio );
    size_t n = adjustable[a];
    CHECK_COMPARE( NoSites[a], <, reservoirs[k][n].size() )( Qn[k].size() )( volume[k] );
    Qn[k].insert( reservoirs[k][n][NoSites[a]++] );
    }
   }
  }
 
 // 3) Generate OneParticleBasis and Subsystems.
 OneParticleBasis QnTot;
 for(size_t k = 0; k < Qn.size(); k++)  QnTot.insert( Qn[k].labels() );
 PRECONDITION( !Nadd_fullSystem || QnTot.size() );
 
 if( !qnLabelNames.size() )  qnLabelNames = qnNames[coord_idcs];
 LabelList<string> qnNames_ll(qnNames);
 Subsystems sys( qnNames.labels(), qnNames_ll(qnLabelNames), QnTot.labels() );
 for(size_t k = 0; k < Qn.size(); k++)
  sys.append_subsystem( Qn[k].labels() );
 
 for(size_t n = 0; n < Nadd_fullSystem; n++)
  sys.append_system();
 
 if( show_subsystems )  { TRACE( sys )( sys.is_complete() ); }
 else  { TRACE( sys.size() )( sys.NoSubsystems() )( sys.is_complete() ); }
 
 
 cout << "\n============================================================== store Subsystems ==\n";
 if( !out_subsys.size() )
  out_subsys = suffix_stripped( pathName_stripped( in_symm ), ".Symm.bin" )
               + "_" + label_subsys + ".Subsys.bin";
 
 store( out_subsys, hist, sys );
 
 
 cout<<endl<<"END"<<endl;
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 }
