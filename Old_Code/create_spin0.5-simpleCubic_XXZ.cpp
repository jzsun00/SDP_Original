#include "algebra/algebra00.h"
#include "io/operator00_io.h"
#include "io/symmetries00_io.h"

#include "lattice/Subsystems00.h"	// just for random operators
#include "latticeAlgebra/operatorTransform00.h"	// just for random operators


#include "fileioL00.h"
#include "linalg/TensorL02.h"


using namespace std;
using namespace LinearAlgebra;


//__________________________________________________________________________________________________
//########################################################################################### MAIN #
int main(int argc, char* argv[])
 {
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 srand(0);
 
 string model   = "spin0.5-simpleCubic_XXZ";
 string program = "create_"+model;
 string Hk      = "H = 0.5*J_{xy}*\\sum_{<ij>}( S^+_i S^-_j + h.c. ) + Jz*\\sum_{<ij>} S^z_i S^z_j - h*\\sum_i S^z_i";
 
 int  algebraF = 2;
 bool realF    = 1;
 
 cout<<" >>> Create spin-1/2 XXZ model on a simple cubic lattice <<< 09.08.2011\n";
 cout<<" usage:  "<< program << endl;
 cout<<"		 -L <PBC system size (0 indicates infinite size)> ...\n";
 cout<<"		(-Jxy <.=1>)  (-Jz <.=1>)  (-h <.=0>)\n";
 cout<<"		(-dont_tell_pointSymmetries)\n";
 cout<<"		(-out_symm     <*.Symm.bin        file =std string>)\n";
 cout<<"		(-out_op       <*.PolyQn.list.bin file =std string>)\n";
 cout<<" Hamiltonian: " << Hk << endl;
 cout<<endl;
 
 sdp_gs::PointI L;
 sdp_gs::Float  Jxy = 1;
 sdp_gs::Float  Jz  = 1;
 sdp_gs::Float  h   = 0;
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
 argID["-Jxy"]		= 10;
 argID["-Jz"]		= 12;
 argID["-h"]		= 14;
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
   case  10:  Jxy  = boost::lexical_cast<sdp_gs::Float>( argv[argn+1] );  argn++;  break;
   case  12:  Jz   = boost::lexical_cast<sdp_gs::Float>( argv[argn+1] );  argn++;  break;
   case  14:  h    = boost::lexical_cast<sdp_gs::Float>( argv[argn+1] );  argn++;  break;
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
 
 string label_latt = "spin0.5-simpleCubic" + sdp_gs::vectorLabel( L, "L" );
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
 PolynomialQn H = -h * (ladder_p(qn0)*ladder_m(qn0)-0.5*Id);
 for(size_t d = 0; d < D; d++)
  H += 0.5*Jxy*( ladder_p(qn0)*ladder_m(qn[d]) + herm(ladder_p(qn0)*ladder_m(qn[d])) )
          + Jz* (ladder_p(qn0)*ladder_m(qn0)-0.5*Id) * (ladder_p(qn[d])*ladder_m(qn[d])-0.5*Id);
 
 TRACE( H );
 H = normalOrdered( H, algebraF );
 CHECK_EQUAL( H, herm(H) );
 
 string param_H = string("XXZ_")
                 + "Jxy" + boost::lexical_cast<string>(Jxy)
                 +  "Jz" + boost::lexical_cast<string>(Jz)
                 +   "h" + boost::lexical_cast<string>(h);
 string label_H = label_latt + "_" + param_H;
 
 
 PolynomialQn total_N, total_Sz;
 if( min(latt.L()) > 0 )	// We consider a finite system.
  {
  cout << "\n=================================================== create conserved quantities ==\n";
  TensorIndex ti(L);
  for(TensorIndexState i = ti.begin(); i.is_valid(); ++i)
   {
   QuantumNumber qni = direct_sum( QuantumNumber(1,1), i() );
   total_N += ladder_p(qni)*ladder_m(qni);
   }
  total_Sz = total_N - (0.5*ti.Dim())*Id;
  CHECK_EQUAL( total_N,  herm(total_N) );
  CHECK_EQUAL( total_Sz, herm(total_Sz) );
  }
 
 
 cout << "\n============================================================= store Hamiltonian ==\n";
 if( !out_op.size() )
  out_op = label_H + ".PolyQn.list.bin";
 
 {
 LabelledList<PolynomialQn> op;
 op["Id"] = Id;
 op["H"]  = H;
 op["Sp"] = ladder_p(qn0);
 op["Sm"] = ladder_m(qn0);
 op["Sz"] = ladder_p(qn0)*ladder_m(qn0)-0.5*Id;
 if( min(latt.L()) > 0 )
  {
  op["tN"]  = total_N;
  op["tNtN"]  = total_N*total_N;
 /* op["tSz"] = total_Sz;
  
  QuantumNumbers qns = sublattice_cube( latt, latt.L() );
  qns = qn_directSum( QuantumNumbers( qns.size(), QuantumNumber(1,1) ), qns );
  Subsystems sys( qnNames.labels(), labelIdcs, qns );
  sys.append_subsystem( qns );
  sys.append_subsystem( qns );
  for(size_t n = 0; n < 10; n++)
   {
   string label = "rand"+boost::lexical_cast<string>(n);
   op[label] = polynomialQn(
                random_Polynomial( sys(), sdp_gs::Indices(sys.NoSubsystems(),latt.dim()) ), sys.basis() );
   op[label] += herm(op[label]);
   }*/
  }
 TRACE( op.labels() );
 
 History hist0 = hist;
 hist0[-1].info_s["label_Lattice"]     = label_latt;
 hist0[-1].info_s["label_Hamiltonian"] = label_H;
 hist0[-1].info_s["param_Hamiltonian"] = param_H;
 hist0[-1].info_s["Hamiltonian_formula"] = Hk;
 PolynomialQn H = -h * (ladder_p(qn0)*ladder_m(qn0)-0.5*Id);
 for(size_t d = 0; d < D; d++)
  H += 0.5*Jxy*( ladder_p(qn0)*ladder_m(qn[d]) + herm(ladder_p(qn0)*ladder_m(qn[d])) )
          + Jz* (ladder_p(qn0)*ladder_m(qn0)-0.5*Id) * (ladder_p(qn[d])*ladder_m(qn[d])-0.5*Id);
 
 hist0[-1].param_d["Jxy"] = Jxy;
 hist0[-1].param_d["Jz"]  = Jz;
 hist0[-1].param_d["h"]   = h;
 
 store( out_op, hist0, op, qnNames );
 }
 
 
 cout<<endl<<"END"<<endl;
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 }
