#include "algebra/algebra00.h"
#include "io/GreenFunction00_io.h"

#include "fileioL00.h"

#include <fstream>
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"


using namespace std;
using namespace LinearAlgebra;



//##################################################################################################

//--------------------------------------------------------------------------[ erase_all ]-
string erase_all_copy(string const& txt, string const& to_remove)
 {
 set<char> rm;
 for(size_t n = 0; n < to_remove.size(); n++)  rm.insert( to_remove[n] );
 string ret( txt );
 size_t cnt = 0;
 for(size_t n = 0; n < txt.size(); n++)
  if( rm.find(txt[n]) == rm.end() )
   ret[cnt++] = txt[n];
 return ret.substr(0,cnt);
 }


//------------------------------------------------------------------------[ read_vector ]-
Vector<sdp_gs::Float> read_vector(string const& txt)
 {
 string clean_txt = erase_all_copy( txt, " \t{()}\n" );
 list<string> split_txt;
 boost::split( split_txt, clean_txt, boost::algorithm::is_any_of(",") );
 
 Vector<sdp_gs::Float> ret( split_txt.size(), 3e8 );
 size_t cnt = 0;
 for(list<string>::const_iterator i = split_txt.begin(); i != split_txt.end(); ++i, ++cnt)
  ret[cnt] = boost::lexical_cast<sdp_gs::Float>( *i );
 
 return ret;
 }



//__________________________________________________________________________________________________
//########################################################################################### MAIN #
int main(int argc, char* argv[])
 {
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 
 string program = "import-G-from-sdpa";
 
 cout<<" >>> Import a Green's function from an SDPA ini file <<< 29.08.2011\n";
 cout<<" usage:  "<< program << endl;
 cout<<"		 -in_ini <*.sdpa.ini or *.sdplr.ini file>  (-file_format <'sdpa' or 'sdplr'>)\n";
 cout<<"		(-in_history <*.history.bin>)  (-in_GreenBasis <*.Gbasis.bin file>)\n";
 cout<<"		(-out_G <*.G.bin =std string>)\n";
 cout<<"  The file format needs only to be specified if the file suffixes are non-standard (and the format is not 'sdpa').\n";
 cout<<endl;
 
 string in_ini;
 string in_ini_format = "sdpa";
 string in_history;
 string in_GreenBasis;
 string out_G;
 
 
 //=================================================================================
 //-[read in command line parameters]-----------------------------------------------
 if( argc<2) return 1;
 
 int argMode = 0;
 set<string> argModable;
 map<string,int> argID;
 argID["-in_ini"]		= 1;		//argModable.insert("-opLabels");
 argID["-in_ini_format"]	= 3;
 argID["-in_history"]		= 5;
 argID["-in_GreenBasis"]	= 7;
 argID["-out_G"]		= 10;
 
 for(int argn=1; argn<argc; argn++)
  {
  int aID;  string argStr( argv[argn] );
  if( argStr[0] == '-' && (argStr[1]<'0'||argStr[1]>'9') )
   { aID = argID[argStr];  argMode = 0;  if(argModable.count(argStr)) argMode = aID; }
  else { aID=argMode;  argn--; }
  switch( aID )
   {
   case   1:  in_ini        = argv[argn+1];  argn++;  break;
   case   3:  in_ini_format = argv[argn+1];  argn++;  break;
   case   5:  in_history    = argv[argn+1];  argn++;  break;
   case   7:  in_GreenBasis = argv[argn+1];  argn++;  break;
   case  10:  out_G         = argv[argn+1];  argn++;  break;
   default :  PANIC("Wrong command line arguments!")(argStr);  break;
   }
  }
 
 PRECONDITION( fileExists(in_ini) )( in_ini );
 if( in_GreenBasis.size() )
  { PRECONDITION( fileExists(in_GreenBasis) )( in_GreenBasis ); }
 if( file_suffixes( in_ini, 2 ) == ".sdplr.ini")  in_ini_format = "sdplr";
 
 
 cout << "\n================================================================== read History ==\n";
 if( !in_history.size() )
  {
  in_history = suffix_stripped( in_ini, 2 ) + ".sdpa.history.bin";
  TRACE( in_history );
  }
 PRECONDITION( fileExists(in_history) )( in_history );
 
 History hist;
 {
 unsigned ID;
 int ver;
 read( in_history, hist, ID, ver);
 CHECK_EQUAL( ID, History_io_ID );
 CHECK_EQUAL( ver, 0 );
 }
 
 int exp_idx = hist.rfind_info( "program", "export-minProb-to-sdpa" );
 CHECK_COMPARE( exp_idx, >=, -1 );
 bool leave_subdir = (bool)hist[exp_idx].get_param_d("create_subdir");
 string file_GreenBasis = find_file( hist[exp_idx].get_param_s("file_GreenBasis"), in_history );
 
 
 cout << "\n=============================================================== read GreenBasis ==\n";
 History    hist_G;
 GreenBasis G_basis;
 LabelList<string> qnNames;
 load_cmp( in_GreenBasis, file_GreenBasis, hist_G, G_basis, qnNames );
 
 
 cout << "\n======================================================================= History ==\n";
 History::Action a;
 a.info_s["program"]      = program;
 a.info_s["command_line"] = commandLine( argc, argv );
 a.info_s["cwd"]          = current_workingPath();
 a.info_s["date"]         = dateTimeString();
 a.param_s["file_sdpa_ini"]        = in_ini;
 a.param_s["file_sdpa_ini_format"] = in_ini_format;
 a.param_s["file_sdpa_history"] = in_history;
 a.param_d["leave_subdir"]      = leave_subdir;
 hist.push_back(a);
 
 
 cout << "\n================================== read Green's function from the sdpa-ini file ==\n";
 ifstream ifs( in_ini.c_str() );
 string temp;
 Vector<sdp_gs::Float> Gvec;
 
 if( in_ini_format == "sdpa" )
  {
  while(1)
   {
   getline( ifs, temp );
   PRECONDITION( !ifs.fail() );
   if( temp.size() > 2 )  break;
   }
  ifs.close();
  CHECK_COMPARE( temp.size(), >, 2 );
  Gvec = read_vector( temp );
  }
 else if( in_ini_format == "sdplr" )
  {
  while(1)
   {
   getline( ifs, temp );
   PRECONDITION( !ifs.fail() );
   if( temp.size() > 2 )  break;
   }
  vector<string> split_txt;
  boost::split( split_txt, temp, boost::algorithm::is_any_of(" ") );
  CHECK_COMPARE( split_txt.size(), >=, 3 );
  CHECK_EQUAL( split_txt[0], "dual" );
  CHECK_EQUAL( split_txt[1], "variable" );
  size_t cnt = boost::lexical_cast<size_t>( split_txt[2] );
  CHECK_EQUAL( cnt+1, G_basis.size() );
  Gvec = Vector<sdp_gs::Float>( cnt );
  for(size_t n = 0; n < cnt; n++)
   {
   getline( ifs, temp );
   PRECONDITION( !ifs.fail() );
   Gvec[n] = -boost::lexical_cast<sdp_gs::Float>( temp );
   }
  }
 else
  { PANIC("Unknown sdp ini file format.")(in_ini_format);  }
 
 GreenFunction G;
 CHECK_EQUAL( G_basis.size(), Gvec.size()+1 );
 CHECK_EQUAL( G_basis.elementsReal()[0], MonomialQn() );
 G.first  = direct_sum( Vector<sdp_gs::Float>(1,1), Gvec[range(0,G_basis.sizeReal()-1)] );
 G.second = Gvec[range(G_basis.sizeReal()-1,Gvec.size())];
 
 
 cout << "\n=========================================================== store GreenFunction ==\n";
 if( !out_G.size() )
  {
  if( leave_subdir )
   {
   size_t pos = in_ini.rfind(".export");
   if( pos == string::npos )  pos = in_ini.rfind("/");
   PRECONDITION( pos != string::npos );
   out_G = in_ini.substr(0,pos) + ".G.bin";
   }
  else
   out_G = suffix_stripped( pathName_stripped( in_ini ), 2 ) + ".G.bin";
  }
 
 store( out_G, hist, G, in_GreenBasis );
 
 
 cout<<endl<<"END"<<endl;
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 }
