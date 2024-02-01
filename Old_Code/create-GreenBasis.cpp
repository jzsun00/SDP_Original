#include "io/operator00_io.h"
#include "io/symmetries00_io.h"
#include "io/GreenFunction00_io.h"

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
 
 rusage usage_start;
 getrusage( RUSAGE_SELF, &usage_start );
 string program = "create-GreenBasis";
 
 cout<<" >>> Create Green's function basis for a certain constraint operator basis, set of Polynomials, and Symmetries <<< 14.08.2011\n";
 cout<<" usage:  "<< program << endl;
 cout<<"		 -in_symm <*.Symm.bin file>  (-override_realF <.=-1>)\n";
 cout<<"		(-in_constrOpBasis <*.PolysQn.bin file>)\n";
 cout<<"		(-set_in_poly <*.PolyQn.list.bin file>  -opLabels <operator labels for the current poly file> ...)\n";
 cout<<"		(-show_history)  (-show_GreenBasis)  (-dont_use_repMap)\n";
 cout<<"		(-label_GreenBasis <label for the resulting Green function basis>)\n";
 cout<<"		(-out_GreenBasis <*.Gbasis.bin file =std string>)\n";
 cout<<" You can call '-set_in_poly' several times. Each time, use '-op_labels' to specify which operators\n";
 cout<<"  from the current file are to be used.\n";
 cout<<" '-override_realF 0': set symm.realF to False, '-override_realF 1': set symm.realF to True.\n"; 
 cout<<endl;
 
 bool show_history    = 0;
 bool show_GreenBasis = 0;
 bool use_repMap      = 1;
 int  override_realF  = -1;
 string in_symm;
 string in_constrOpBasis;
 Vector<string>          in_poly;
 Vector<Vector<string> > in_opLabel;
 string out_GreenBasis;
 string label_GreenBasis;
 
 
 //=================================================================================
 //-[read in command line parameters]-----------------------------------------------
 if( argc<2) return 1;
 
 int argMode = 0;
 set<string> argModable;
 map<string,int> argID;
 argID["-in_symm"]		= 1;
 argID["-override_realF"]	= 3;
 argID["-in_constrOpBasis"]	= 10;
 argID["-set_in_poly"]		= 20;
 argID["-opLabels"]		= 22;		argModable.insert("-opLabels");
 argID["-show_history"]		= 30;
 argID["-show_GreenBasis"]	= 34;
 argID["-dont_use_repMap"]	= 36;
 argID["-label_GreenBasis"]	= 40;
 argID["-out_GreenBasis"]	= 50;
 
 for(int argn=1; argn<argc; argn++)
  {
  int aID;  string argStr( argv[argn] );
  if( argStr[0] == '-' && (argStr[1]<'0'||argStr[1]>'9') )
   { aID = argID[argStr];  argMode = 0;  if(argModable.count(argStr)) argMode = aID; }
  else { aID=argMode;  argn--; }
  switch( aID )
   {
   case   1:  in_symm = argv[argn+1];  argn++;  break;
   case   3:  override_realF = boost::lexical_cast<int>( argv[argn+1] );  argn++;  break;
   case  10:  in_constrOpBasis = argv[argn+1];  argn++;  break;
   case  20:  in_poly << argv[argn+1];  in_opLabel << Vector<string>();  argn++;  break;
   case  22:  in_opLabel[ in_opLabel.size()-1 ] << argv[argn+1];  argn++;  break;
   case  30:  show_history    = 1;  break;
   case  34:  show_GreenBasis = 1;  break;
   case  36:  use_repMap      = 0;  break;
   case  40:  label_GreenBasis = argv[argn+1];  argn++;  break;
   case  50:  out_GreenBasis   = argv[argn+1];  argn++;  break;
   default :  PANIC("Wrong command line arguments!")(argStr);  break;
   }
  }
 
 RANGE_CHECK( override_realF, -1, 1 );
 PRECONDITION( fileExists(in_symm) )( in_symm );
 PRECONDITION( in_poly.size() || in_constrOpBasis.size() );
 PRECONDITION( !in_constrOpBasis.size() || fileExists(in_constrOpBasis) )( in_constrOpBasis );
 for(size_t n = 0; n < in_poly.size(); n++)
  { PRECONDITION( fileExists( in_poly[n] ) )( in_poly[n] ); }
 
 
 cout << "\n========================= read Symmetries, Polynomials, and contraint operators ==\n";
 // Symmetries
 History    hist_symm;
 Symmetries symm;
 LabelList<string> qnNames;
 read( in_symm, hist_symm, symm, qnNames );
 int algebraF = (int)hist_symm[-1].get_info_d("algebraF");
 
 if( show_history )  { TRACE( hist_symm ); }
 
 // constraint operator basis
 History hist_constrOpBasis;
 PolynomialsQn constrOpBasis;
 if( in_constrOpBasis.size() )
  {
  LabelList<string> qnNames2;
  read( in_constrOpBasis, hist_constrOpBasis, constrOpBasis, qnNames2 );
  if( show_history )  { TRACE( hist_constrOpBasis ); }
  CHECK_EQUAL( qnNames2, qnNames )( in_constrOpBasis );
  CHECK_EQUAL( (int)hist_constrOpBasis[-1].get_info_d("algebraF"), algebraF );
  }
 
 // Polynomials
 History       hist_poly;
 PolynomialsQn op;
 for(size_t n = 0; n < in_poly.size(); n++)
  {
  LabelledList<PolynomialQn> list;
  LabelList<string> qnNames2;
  read( in_poly[n], hist_poly, list, qnNames2 );
  if( show_history )  { TRACE( hist_poly ); }
  CHECK_EQUAL( qnNames2, qnNames )( in_poly[n] );
  CHECK_EQUAL( (int)hist_poly[-1].get_info_d("algebraF"), algebraF );
  if( override_realF == -1 )
   { CHECK_EQUAL( (int)hist_poly[-1].get_info_d("realF"), symm.is_real() ); }
  
  for(size_t nOp = 0; nOp < in_opLabel[n].size(); nOp++)
   {
   PRECONDITION( list.contains( in_opLabel[n][nOp] ) );
   op << list.find( in_opLabel[n][nOp] )->second;
   }
  }
 
 
 cout << "\n======================================================================= History ==\n";
 History hist;
 if( in_constrOpBasis.size() )  hist = hist_constrOpBasis;
 else  hist = hist_poly;
 
 History::Action a;
 a.info_s["program"]      = program;
 a.info_s["command_line"] = commandLine( argc, argv );
 a.info_s["cwd"]          = current_workingPath();
 a.info_s["date"]         = dateTimeString();
 a.param_s["file_Symmetries"]    = in_symm;
 a.param_s["file_constrOpBasis"] = in_constrOpBasis;
 hist.push_back(a);
 
 if( label_GreenBasis.size() )
  hist[-1].info_s["label_GreenBasis"] = label_GreenBasis;
 
 
 cout << "\n=================================================== create Green function basis ==\n";
 if( override_realF != -1 )
  symm.set_realF( bool(override_realF) );
 
 hist[-1].param_d["override_realF"] = override_realF;
 hist[-1].param_d["realF"]          = symm.is_real();
 
 GreenBasis G_basis( algebraF, OneParticleBasis(), symm );
 
 if( constrOpBasis.size() )
  {
  cout << "\n Collecting elements needed for positivity constraints ...\n";
  TRACE( constrOpBasis.size() );
  size_t m0 = 1;
  for(size_t m = 0; m < constrOpBasis.size(); m++)
   {
   if( m == m0 )  { cout << " " << m << flush;  m0*=2; }
   for(size_t n = 0; n < constrOpBasis.size(); n++)
    G_basis.collect_elements( herm(constrOpBasis[m])*constrOpBasis[n], use_repMap );
   }
  
  TRACE( G_basis.sizeReal() )( G_basis.sizeImag() )( G_basis.representativeMap().size() )
       ( G_basis.distance_min_Re() )( G_basis.distance_min_Im() );
  }
 
 if( op.size() )
  {
  cout << "\n Collecting elements needed for specified operators ...\n";
  for(size_t n = 0; n < op.size(); n++)
   G_basis.collect_elements( op[n], use_repMap );
  
  TRACE( G_basis.sizeReal() )( G_basis.sizeImag() )( G_basis.representativeMap().size() )
       ( G_basis.distance_min_Re() )( G_basis.distance_min_Im() );
  }
 
 if( show_GreenBasis )  { TRACE( G_basis ); }
 
 
 cout << "\n==================================================== store Green function basis ==\n";
 if( !out_GreenBasis.size() && in_constrOpBasis.size() )
  out_GreenBasis = suffix_stripped( pathName_stripped( in_constrOpBasis ), ".PolysQn.bin" )
                   + (label_GreenBasis.size() ? ("_" + label_GreenBasis) : "") + ".Gbasis.bin";
 if( !out_GreenBasis.size() && in_poly.size() )
  out_GreenBasis = suffix_stripped( pathName_stripped( in_poly[0] ), ".PolyQn.list.bin" )
                   + (label_GreenBasis.size() ? ("_" + label_GreenBasis) : "") + ".Gbasis.bin";
 PRECONDITION( out_GreenBasis.size() );
 
 hist[-1].info_d["time_cpu_used"] = cpuTime_used( usage_start );
 store( out_GreenBasis, hist, G_basis, qnNames );
 
 
 cout<<endl<<"END"<<endl;
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 }
