#include "latticeAlgebra/operatorBasis00.h"
#include "latticeAlgebra/operatorTransform00.h"
#include "io/Subsystems00_io.h"
#include "io/operator00_io.h"

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
 
 string program = "create-opBasis";
 
 cout<<" >>> Create operator basis in form of a PolynomialsQn list <<< 13.08.2011\n";
 cout<<" usage:  "<< program << endl;
 cout<<"		 -create_full <*.Subsys.bin file>  (-show_history)  (-show_opBasis)  (-show_opBasisQn)\n";
 cout<<"		(-label_opBasis <label for the resulting operator basis>)\n";
 cout<<"		(-out_opBasis <*.PolysQn.bin file =std string>)\n";
 cout<<endl;
 
 bool do_full = 0;
 bool show_history   = 0;
 bool show_opBasis   = 0;
 bool show_opBasisQn = 0;
 string in_subsys;
 string out_opBasis;
 string label_opBasis;
 
 
 //=================================================================================
 //-[read in command line parameters]-----------------------------------------------
 if( argc<2) return 1;
 
 int argMode = 0;
 set<string> argModable;
 map<string,int> argID;
 argID["-create_full"]		= 1;
 argID["-show_history"]		= 10;
 argID["-show_opBasis"]		= 12;
 argID["-show_opBasisQn"]	= 14;
 argID["-label_opBasis"]	= 20;		//argModable.insert("-extra_param");
 argID["-out_basis"]		= 30;
 
 for(int argn=1; argn<argc; argn++)
  {
  int aID;  string argStr( argv[argn] );
  if( argStr[0] == '-' && (argStr[1]<'0'||argStr[1]>'9') )
   { aID = argID[argStr];  argMode = 0;  if(argModable.count(argStr)) argMode = aID; }
  else { aID=argMode;  argn--; }
  switch( aID )
   {
   case   1:  do_full = 1;
              in_subsys = argv[argn+1];  argn++;  break;
   case  10:  show_history = 1;  break;
   case  12:  show_opBasis   = 1;  break;
   case  14:  show_opBasisQn = 1;  break;
   case  20:  label_opBasis = argv[argn+1];  argn++;  break;
   case  30:  out_opBasis   = argv[argn+1];  argn++;  break;
   default :  PANIC("Wrong command line arguments!")(argStr);  break;
   }
  }
 
 PRECONDITION( do_full );	// TODO: Add some alternatives later?
 PRECONDITION( fileExists(in_subsys) )( in_subsys );
 
 
 cout << "\n========================================================== read Subsystems file ==\n";
 History hist;
 Subsystems sys;
 read( in_subsys, hist, sys );
 
 int algebraF = (int)hist[-1].get_info_d("algebraF");
 string label_subsys = hist[-1].get_info_s("label_Subsystems");
 
 if( show_history )  { TRACE( hist ); }
 
 
 cout << "\n======================================================================= History ==\n";
 History::Action a;
 a.info_s["program"]      = program;
 a.info_s["command_line"] = commandLine( argc, argv );
 a.info_s["cwd"]          = current_workingPath();
 a.info_s["date"]         = dateTimeString();
 a.info_d["algebraF"]     = algebraF;
 a.param_s["file_Subsystems"] = in_subsys;
 hist.push_back(a);
 
 if( !label_opBasis.size() )  label_opBasis = "full";
 hist[-1].info_s["label_constraintOperatorBasis"] = label_opBasis;
 
 
 cout << "\n========================================================= create operator basis ==\n";
 sdp_gs::Indices NoElements;
 Monomials opBasis = operatorBasis_full( sys, algebraF, NoElements );
 
 TRACE( algebraF )( NoElements );
 if( show_opBasis )    { TRACE( opBasis ); }
 if( show_opBasisQn )  { TRACE( monomialsQn( opBasis, sys.basis() ) ); }
 
 
 cout << "\n========================================================== store operator basis ==\n";
 if( !out_opBasis.size() )
  out_opBasis = suffix_stripped( pathName_stripped( in_subsys ), ".Subsys.bin" )
                + "_" + label_opBasis + ".PolysQn.bin";
 
 store( out_opBasis, hist, polynomialsQn( opBasis, sys.basis() ), sys.qnNames() );
 
 
 cout<<endl<<"END"<<endl;
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 }
