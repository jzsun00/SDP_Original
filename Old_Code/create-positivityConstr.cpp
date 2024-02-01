#include "io/operator00_io.h"
#include "io/GreenFunction00_io.h"
#include "io/PositivityConstr01_io.h"
#include "constraints/sdpInfo00.h"

#include "fileioL00.h"
#include "parseMathL03.h"


using namespace std;
using namespace LinearAlgebra;



//##################################################################################################

//--------------------------------------------------------------[ do_randSpec_and_store ]-
template <typename T>
void do_randSpec_and_store(T const& Mb, Vector<sdp_gs::Indices> const& blocks_global,
                           bool show_random_spectra, bool probe_hermiticity,
                           GreenBasis const& G_basis,
                           string const& out_positivityConstr, History const& hist,
                           string const& in_opBasis, string const& in_GreenBasis)
 {
 using namespace LinearAlgebra;
 
 if( probe_hermiticity )
  {
  cout << "\n========================================= probing constraint matrix hermiticity ==\n";
  check_constraintHermiticity( Mb );
  }
 
 if( show_random_spectra )
  {
  cout << "\n============================================== show random G constraint spectra ==\n";
  show_constraintRandomSpectra( cout, G_basis, Mb );
  }
 
 cout << "\n================================================== store positivity constraints ==\n";
 store( out_positivityConstr, hist, Mb, blocks_global, in_opBasis, in_GreenBasis );
 }



//__________________________________________________________________________________________________
//########################################################################################### MAIN #
int main(int argc, char* argv[])
 {
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 
 rusage usage_start;
 getrusage( RUSAGE_SELF, &usage_start );
 string program = "create-positivityConstr";
 
 cout<<" >>> Create positivity constraints for a certain constraint operator basis, and Green basis <<< 15.08.2011\n";
 cout<<" usage:  "<< program << endl;
 cout<<"		 -in_GreenBasis <*.Gbasis.bin file>  -in_opBasis <*.PolysQn.bin file>  (-show_history)\n";
 cout<<"                (-dont_make_real)  (-dont_separate_blocks)  (-show_blockElements)  (-show_blockElementReps <qnIdx> ...)  (-show_nnz)\n";
 cout<<"                (-force_separation <formula for qn condition> ...)\n";
 cout<<"		(-show_random_spectra)  (-dont_probe_hermiticity)\n";
 cout<<"		(-label_positivityConstr <label for the resulting constraints>)\n";
 cout<<"		(-out_positivityConstr <*.Constr.bin file =std string>)\n";
 cout<<" The formula for the quantum number condition may depend on all quantum number names,\n";
 cout<<"  and 'degree' and 'degree_min', the polynomial degrees.\n";
 cout<<endl;
 
 bool show_history = 0;
 bool make_real    = 1;
 bool do_separate_blocks    = 1;
 bool show_blockElements    = 0;
 bool show_blockElementReps = 0;
 bool show_nnz     = 0;
 bool show_random_spectra = 0;
 bool probe_hermiticity   = 1;
 string in_GreenBasis;
 string in_opBasis;
 string out_positivityConstr;
 string label_positivityConstr;
 sdp_gs::Indices qn_idcs;
 Vector<ParseMath::Formula> sepCond;
 
 
 //=================================================================================
 //-[read in command line parameters]-----------------------------------------------
 if( argc<2) return 1;
 
 int argMode = 0;
 set<string> argModable;
 map<string,int> argID;
 argID["-in_GreenBasis"]	= 1;		//argModable.insert("-opLabels");
 argID["-in_opBasis"]		= 3;
 argID["-show_history"]		= 20;
 argID["-dont_make_real"]	= 30;
 argID["-dont_separate_blocks"]	= 32;
 argID["-force_separation"]	= 33;		argModable.insert("-force_separation");
 argID["-show_blockElements"]	= 34;
 argID["-show_blockElementReps"]   = 35;	argModable.insert("-show_blockElementReps");
 argID["-show_nnz"]		   = 36;
 argID["-show_random_spectra"]     = 38;
 argID["-dont_probe_hermiticity"]  = 40;
 argID["-label_positivityConstr"]  = 50;
 argID["-out_positivityConstr"]	   = 60;
 
 for(int argn=1; argn<argc; argn++)
  {
  int aID;  string argStr( argv[argn] );
  if( argStr[0] == '-' && (argStr[1]<'0'||argStr[1]>'9') )
   { aID = argID[argStr];  argMode = 0;  if(argModable.count(argStr)) argMode = aID; }
  else { aID=argMode;  argn--; }
  switch( aID )
   {
   case   1:  in_GreenBasis = argv[argn+1];  argn++;  break;
   case   3:  in_opBasis    = argv[argn+1];  argn++;  break;
   case  20:  show_history = 1;  break;
   case  30:  make_real    = 0;  break;
   case  32:  do_separate_blocks  = 0;  break;
   case  33:  sepCond << ParseMath::parse( argv[argn+1] );  argn++;  break;
   case  34:  show_blockElements  = 1;  break;
   case  35:  show_blockElementReps = 1;  PRECONDITION( argn+1 < argc );
              qn_idcs << boost::lexical_cast<size_t>( argv[argn+1] );  argn++;  break;
   case  36:  show_nnz     = 1;  break;
   case  38:  show_random_spectra = 1;  break;
   case  40:  probe_hermiticity   = 0;  break;
   case  50:  label_positivityConstr = argv[argn+1];  argn++;  break;
   case  60:  out_positivityConstr   = argv[argn+1];  argn++;  break;
   default :  PANIC("Wrong command line arguments!")(argStr);  break;
   }
  }
 
 PRECONDITION( fileExists(in_GreenBasis) )( in_GreenBasis );
 PRECONDITION( fileExists(in_opBasis) )   ( in_opBasis );
 
 
 cout << "\n================================ read GreenBasis, and constraint operator basis ==\n";
 // GreenBasis
 History    hist_G;
 GreenBasis G_basis;
 LabelList<string> qnNames;
 read( in_GreenBasis, hist_G, G_basis, qnNames );
 sdp_gs::Indices qn_idcs0( range(0,qnNames.size()) );
 
 if( show_history )  { TRACE( hist_G ); }
 
 // constraint operator basis
 History hist_opBasis;
 PolynomialsQn opBasis;
 {
 LabelList<string> qnNames2;
 read( in_opBasis, hist_opBasis, opBasis, qnNames2 );
 if( show_history )  { TRACE( hist_opBasis ); }
 CHECK_EQUAL( qnNames2, qnNames );
 CHECK_EQUAL( (int)hist_opBasis[-1].get_info_d("algebraF"), G_basis.algebraF() );
 }
 
 
 cout << "\n======================================================================= History ==\n";
 History hist( hist_opBasis );
 hist.push_back( hist_G[-1] );
 
 History::Action a;
 a.info_s["program"]      = program;
 a.info_s["command_line"] = commandLine( argc, argv );
 a.info_s["cwd"]          = current_workingPath();
 a.info_s["date"]         = dateTimeString();
 a.param_d["make_real"]       = make_real;
 a.param_d["separate_blocks"] = do_separate_blocks;
 hist.push_back(a);
 
 if( !label_positivityConstr.size() )
  label_positivityConstr =   "re"  + boost::lexical_cast<string>(make_real)
                           + "sep" + boost::lexical_cast<string>(do_separate_blocks);
 hist[-1].info_s["label_PositivityConstraints"] = label_positivityConstr;
 
 if( !out_positivityConstr.size() )
  out_positivityConstr = suffix_stripped( pathName_stripped( in_opBasis ), ".PolysQn.bin" )
                         + "_" + label_positivityConstr + ".Constr.bin";
 
 
 cout << "\n================================================= create positivity constraints ==\n";
 /*Polynomials OpBasis;
 if( show_blockElements || show_blockElementReps || sepCond.size() )
  {
  OneParticleBasis basis;
  collect_oneParticleBasis( basis, opBasis );
  OpBasis = polynomials( opBasis, basis );
  }*/
 
 cout << "\n------------- create full positivity constraints --\n";
 Vector<PositivityConstraintsCF> Mca( 1 );
 Mca[0] = PositivityConstraintsCF::from_rootBasis( opBasis, G_basis );
 cout << endl;
 
 CHECK_COMPARE( Mca[0].normInf_imag_Re(), <, 100*zeroThresh );
 CHECK_COMPARE( Mca[0].normInf_real_Im(), <, 100*zeroThresh );
 
 Vector<sdp_gs::Indices> blocks_a;
 if( !sepCond.size() )  blocks_a << range(0,Mca[0].dim());
 else
  {
  cout << "\n------------- force user specified block separations --\n";
  Vector<sdp_gs::Indices> blocks_a0( sepCond.size()+1 );
  for(size_t nOp = 0; nOp < opBasis.size(); nOp++)
   {
   ParseMath::VariableMap varMap;
   varMap["degree"]     = degree( opBasis[nOp] );
   varMap["degree_min"] = degree_min( opBasis[nOp] );
   QuantumNumber qn = quantumNumber( opBasis[nOp], qn_idcs0 );
   if( qn.size() )
    for(size_t k = 0; k < qn.size(); k++)  varMap[ qnNames[k] ] = qn[k];
   int nCond0 = -1;
   for(size_t nCond = 0; nCond < sepCond.size(); nCond++)
    {
    ParseMath::Formula cond = ParseMath::eval( sepCond[nCond], varMap );
    if( ParseMath::is_double( cond ) && boost::get<double>(cond) )
     {
     CHECK_EQUAL( nCond0, -1 )( nCond )( sepCond )( opBasis[nOp] );
     nCond0 = nCond;
     }
    }
   blocks_a0[ nCond0+1 ] << nOp;
   }
  for(size_t b = 0; b < blocks_a0.size(); b++)
   if( blocks_a0.size() )  blocks_a << blocks_a0[b];
  if( blocks_a != Vector<sdp_gs::Indices>( 1, range(0,Mca[0].dim()) ) )
   Mca = separate_blocks( Mca[0], blocks_a );
  }
 
 cout << "\n------------- find block structure --\n";
 Vector<Vector<sdp_gs::Indices> > blocks( blocks_a.size() );
 Vector<sdp_gs::Indices> blocks_global;
 for(size_t B = 0; B < blocks.size(); B++)
  {
  if( do_separate_blocks )  blocks[B] = blockStructure( Mca[B] );
  else                      blocks[B] << range(0,Mca[B].dim());
  for(size_t b = 0; b < blocks[B].size(); b++)
   blocks_global << blocks_a[B][ blocks[B][b] ];
  }
 
 show_constraintBlocks( cout << endl, blocks_global );
 if( show_blockElementReps )
  show_constraintBlocks_elementReps( cout << endl, blocks_global, opBasis, qn_idcs );
 if( show_blockElements )
  show_constraintBlocks_elements(    cout << endl, blocks_global, opBasis );
 
 cout << "\n------------- separate blocks --\n";
 Vector<PositivityConstraintsCF> Mcb( blocks_global.size() );
 size_t nb = 0;
 for(size_t B = 0; B < blocks.size(); B++)
  {
  Mcb[ range( nb, nb+blocks[B].size() ) ] = separate_blocks( Mca[B], blocks[B] );
  nb += blocks[B].size();
  Mca[B] = PositivityConstraintsCF();
  }
 
 if( show_nnz ) 
  show_constraintBlocks( cout << endl, Mcb, show_nnz );
 
 hist[-1].info_d["time_cpu_used"] = cpuTime_used( usage_start );
 
 bool do_realify = ( make_real && sdp_gs::my_i != sdp_gs::CFloat(0) );
 if( do_realify )
  {
  cout << "\n------------- making everything real --\n";
  Vector<PositivityConstraintsF> Mb( Mcb.size() );
  for(size_t b = 0; b < Mb.size(); b++)
   {
   Mb[b] = real( Mcb[b] );
   Mcb[b] = PositivityConstraintsCF();	// free memory
   }
  show_constraintBlocks( cout << endl, Mb );
  do_randSpec_and_store( Mb,  blocks_global, show_random_spectra, probe_hermiticity, G_basis,
                         out_positivityConstr, hist, in_opBasis, in_GreenBasis );
  }
 else
  {
  do_randSpec_and_store( Mcb, blocks_global, show_random_spectra, probe_hermiticity, G_basis,
                         out_positivityConstr, hist, in_opBasis, in_GreenBasis );
  }
 
 
 cout<<endl<<"END"<<endl;
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 }
