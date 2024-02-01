#include "algebra/algebra00.h"
#include "io/operator00_io.h"
#include "io/symmetries00_io.h"

#include "fileioL00.h"


using namespace std;
using namespace LinearAlgebra;


//__________________________________________________________________________________________________
//########################################################################################### MAIN #
int main(int argc, char* argv[])
 {
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 
 string model   = "fermi-simpleCubic_Hubbard";
 string program = "create_"+model;
 string Hk      = "H = -t*\\sum_{<ij>}( c^+_i c^-_j + h.c. ) + U*\\sum_{<ij>}(n_i-0.5)*(n_j-0.5) - mu*\\sum_i n_i";
 
 int  algebraF = 1;
 bool realF    = 1;
 
 cout<<" >>> Create spinless Fermi-Hubbard model on a simple cubic lattice <<< 19.08.2011\n";
 cout<<" usage:  "<< program << endl;
 cout<<"		 -L <PBC system size (0 indicates infinite size)> ...\n";
 cout<<"		(-t <.=1>)  (-U <.=1>)  (-mu <.=0>)\n";
 cout<<"		(-dont_tell_pointSymmetries)\n";
 cout<<"		(-out_symm     <*.Symm.bin        file =std string>)\n";
 cout<<"		(-out_op       <*.PolyQn.list.bin file =std string>)\n";
 cout<<" Hamiltonian: " << Hk << endl;
 cout<<endl;
 
 sdp_gs::PointI L;
 sdp_gs::Float  t  = 1;
 sdp_gs::Float  U  = 1;
 sdp_gs::Float  mu = 0;
 bool           tell_pointSymmetries = 1;
 string		out_symm;
 string		out_op;
 
 
 //=================================================================================
 //-[read in command line parameters]-----------------------------------------------
 if( argc<2) return 1;
 
 int argMode = 0;
 set<string> argModable;
 map<string,int> argID;
 argID["-L"]		= 1;		argModable.insert("-L");
 argID["-t"]		= 10;
 argID["-U"]		= 12;
 argID["-mu"]		= 14;
 argID["-dont_tell_pointSymmetries"]	= 20;
 argID["-out_symm"]	= 40;
 argID["-out_op"]	= 42;
 
 for(int argn=1; argn<argc; argn++)
  {
  int aID;  string argStr( argv[argn] );
  if( argStr[0] == '-' && (argStr[1]<'0'||argStr[1]>'9') )
   { aID = argID[argStr];  argMode = 0;  if(argModable.count(argStr)) argMode = aID; }
  else { aID=argMode;  argn--; }
  switch( aID )
   {
   case   1:  L << boost::lexical_cast<int>( argv[argn+1] );  argn++;  break;
   case  10:  t  = boost::lexical_cast<sdp_gs::Float>( argv[argn+1] );  argn++;  break;
   case  12:  U  = boost::lexical_cast<sdp_gs::Float>( argv[argn+1] );  argn++;  break;
   case  14:  mu = boost::lexical_cast<sdp_gs::Float>( argv[argn+1] );  argn++;  break;
   case  20:  tell_pointSymmetries = 0;  break;
   case  40:  out_symm =  argv[argn+1];  argn++;  break;
   case  42:  out_op   =  argv[argn+1];  argn++;  break;
   default :  PANIC("Wrong command line arguments!")(argStr);  break;
   }
  }
 
 PRECONDITION( L.size() );
 
 
 cout << "\n======================================================================= History ==\n";
 History hist;
 History::Action a;
 a.info_s["program"]      = program;
 a.info_s["command_line"] = commandLine( argc, argv );
 a.info_s["cwd"]          = current_workingPath();
 a.info_s["date"]         = dateTimeString();
 a.info_d["algebraF"]     = algebraF;
 a.info_d["realF"]        = realF;
 a.param_d["tell_pointSymmetries"] = tell_pointSymmetries;
 hist.push_back(a);
 
 
 cout << "\n==================================================== create Lattice and qnNames ==\n";
 Lattice latt(L);
 size_t D = L.size();
 
 LabelList<string> qnNames;
 qnNames.push_back("N");
 for(size_t d = 0; d < D; d++)  qnNames.push_back( sdp_gs::dimLabel[d] );
 sdp_gs::Indices labelIdcs(range(1,D+1));

 TRACE( latt )( qnNames )( qnNames[labelIdcs] );
 
 string label_latt = "fermi-simpleCubic" + sdp_gs::vectorLabel( L, "L" );
 for(size_t i = 0; i < L.size(); i++)
  hist[-1].param_d["L"+sdp_gs::dimLabel[i]] = L[i];
 
 
 cout << "\n============================================================= create Symmetries ==\n";
 typedef Symmetries::groupRep_type groupRep_type;
 
 groupRep_type pointGroup = symmGroup_point_cubic( D );
 //TRACE( pointGroup.order() )( pointGroup );
 
 Symmetries symm( realF, Symmetries::lattice_type( range(1,D+1), latt ),
                  Symmetries::groups_type(),
                  sdp_gs::Indices(1,0) );
 Symmetries symm0( symm );
 if( tell_pointSymmetries )
  symm.append_group( IndexSymmetry<groupRep_type>( sdp_gs::Indices(range(1,D+1)), pointGroup ) );
 //symm.append_group( IndexSymmetry<groupRep_type>( sdp_gs::Indices(1,D+1),        spinFlipGroup ) );
 
 TRACE( symm );
 
 
 cout << "\n============================================================== store Symmetries ==\n";
 if( !out_symm.size() )
  out_symm = label_latt + (tell_pointSymmetries?"":"_silent") + ".Symm.bin";
 
 {
 History hist0 = hist;
 hist0[-1].info_s["label_Lattice"] = label_latt;
 store( out_symm, hist0, symm, qnNames );
 }
 
 
 cout << "\n============================================================ create Hamiltonian ==\n";
 QuantumNumber qn0 = direct_sum( QuantumNumber(1,1), QuantumNumber(D,0) );
 QuantumNumbers qn(D,qn0);
 for(size_t d = 0; d < D; d++)  qn[d][d+1] = 1;
 
 MonomialQn Id;
 PolynomialQn H = -mu * (ladder_p(qn0)*ladder_m(qn0));
 for(size_t d = 0; d < D; d++)
  H += -t*( ladder_p(qn0)*ladder_m(qn[d]) + herm(ladder_p(qn0)*ladder_m(qn[d])) )
      + U * (ladder_p(qn0)*ladder_m(qn0)-0.5*Id) * (ladder_p(qn[d])*ladder_m(qn[d])-0.5*Id);
 
 TRACE( H );
 H = normalOrdered( H, algebraF );
 CHECK_EQUAL( H, herm(H) );
 
 string param_H = string("Hubbard_")
                 + "t"  + boost::lexical_cast<string>(t)
                 + "U"  + boost::lexical_cast<string>(U)
                 + "mu" + boost::lexical_cast<string>(mu);
 string label_H = label_latt + "_" + param_H;
 
 
 cout << "\n============================================================= store Hamiltonian ==\n";
 if( !out_op.size() )
  out_op = label_H + ".PolyQn.list.bin";
 
 {
 LabelledList<PolynomialQn> op;
 op["Id"] = Id;
 op["H"]  = H;
 op["ap"] = ladder_p(qn0);
 op["am"] = ladder_m(qn0);
 op["n"]  = ladder_p(qn0)*ladder_m(qn0);
 op["Sz"] = op["n"]-0.5*Id;
 TRACE( op.labels() );
 
 History hist0 = hist;
 hist0[-1].info_s["label_Lattice"]     = label_latt;
 hist0[-1].info_s["label_Hamiltonian"] = label_H;
 hist0[-1].info_s["param_Hamiltonian"] = param_H;
 hist0[-1].info_s["Hamiltonian_formula"] = Hk;
 hist0[-1].param_d["t"]  = t;
 hist0[-1].param_d["U"]  = U;
 hist0[-1].param_d["mu"] = mu;
 
 store( out_op, hist0, op, qnNames );
 }
 
 
 cout<<endl<<"END"<<endl;
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 }
