#include "algebra/algebra00.h"
#include "io/operator00_io.h"
#include "io/GreenFunction00_io.h"
#include "io/PositivityConstr01_io.h"

#include "fileioL00.h"
#include "linearalgebra/eigen.h"	// only needed for "-show_random_spectra"

#include <fstream>
#include "boost/algorithm/string/replace.hpp"


using namespace std;
using namespace LinearAlgebra;



//##################################################################################################

//-------------------------------------------------------------------------[ wrap_lines ]-
string wrap_lines(string const& txt, size_t lineLen_max)
 {
 size_t n0=0;
 string ret;
 while( n0 < txt.size() )
  {
  size_t nn = txt.find( "\n", n0 );
  if( nn == string::npos )  nn = txt.size();
  //TRACE(nn)( nn - n0 )( lineLen_max );
  if( nn - n0 > lineLen_max )
   {
   size_t n_end = std::min( txt.size(), n0 + lineLen_max );
   size_t nr_s = txt.rfind( " ",  n_end );
   size_t nr_t = txt.rfind( "\t", n_end );
   if( nr_s == string::npos )  nr_s = 0;
   if( nr_t == string::npos )  nr_t = 0;
   nr_s = std::max( nr_s, nr_t );
   size_t n1 = nr_s+1;
   if( nr_s <= n0 )  { nr_s = n_end;  n1 = nr_s; }
   //TRACE(n0)( txt.substr( n0, nr_s-n0 ) );
   ret += txt.substr( n0, nr_s-n0 ) + "\\\n";
   n0 = n1;
   }
  else
   {
   //TRACE(n0)( txt.substr( n0, nn+1-n0 ) );
   ret += txt.substr( n0, nn+1-n0 );
   n0 = nn+1;
   }
  }
 return ret;
 }


//-----------------------------------------------------------------------[ print_matrix ]-
// string print_matrix(bool make_denseF, size_t nv, size_t nb,
//                     PositivityConstraintsF::matrix_type const& M)
//  {
//  typedef PositivityConstraintsF::matrix_type matrix_type;
//  ostringstream oss;
//  if( make_denseF )
//   {
//   oss << sdp_gs::MatrixF( M );
//   return oss.str();
//   }
//  for(const_iterator<matrix_type>::type it = iterate(M); it; ++it)
//   for(const_inner_iterator<matrix_type>::type itIn = iterate(it); itIn; ++itIn)
//    if( itIn.index1() <= itIn.index2() )
//     oss << nv << "\t" << nb+1 << "\t"
//         << itIn.index1()+1 << "\t" << itIn.index2()+1 << "\t" << *itIn << endl;
//  string ret = oss.str();
//  if( ret.size() )  ret = ret.substr( 0, ret.size()-1 );
//  return ret;
//  }

string print_matrix(bool make_denseF, size_t nv, size_t nb,
                    PositivityConstraintsF::matrix_type const& M)
 {
 typedef PositivityConstraintsF::matrix_type matrix_type;
 ostringstream oss;
 if( make_denseF )
  {
  oss << make_dense( M );
  return oss.str();
  }
 for(matrix_type::const_iterator i = M.begin(); i != M.end(); ++i)
  if( i->first.first <= i->first.second )
    oss << nv << "\t" << nb+1 << "\t"
        << i->first.first+1 << "\t" << i->first.second+1 << "\t" << i->second << endl;
 string ret = oss.str();
 if( ret.size() )  ret = ret.substr( 0, ret.size()-1 );
 return ret;
 }



//__________________________________________________________________________________________________
//########################################################################################### MAIN #
int main(int argc, char* argv[])
 {
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 
 string program = "export-minProb-to-sdpa";
 
 cout<<" >>> Export a minimization problem with positivity constraints to the file format of SDPA <<< 16.08.2011\n";
 cout<<" usage:  "<< program << endl;
 cout<<"		 -in_positivityConstr <*.Constr.bin file>  (-in_GreenBasis <*.Gbasis.bin file>)\n";
 cout<<"		(-set_in_poly <*.PolyQn.list.bin file>)  (-show_history)\n";
 cout<<"		 -costOp <operator label> \n";
 cout<<"                (-constrain_eigen  <operator label> <value> ...)  (-label_eigen <label for the eigen value constraints>)\n";
 cout<<"                (-constrain_expect <operator label> <value> ...)  (-label_expect <label for the expectation value constraints>)\n";
 cout<<"                (-make_dense)  (-export_history)  (-prec <number of digits =12>)  (-prec_scientific)  (-prec_fixed)\n";
 cout<<"		(-out_sdpa <*.sdpa.dat or *.sdpa.dat-s file =std string>)  (-dont_create_history_file)\n";
 cout<<"		(-create_subdir)  (-create_subdirSuffixed <suffix to add to an export path>)\n";
 cout<<" The given positivity constraints must be real (already), i.e., the *.Constr.bin file should contain\n";
 cout<<"  a Vector<PositivityConstraintsF>.\n";
 cout<<endl;
 
 typedef pair<string,string> OpSpec;	// <PolyQn file name, operator label>
 
 bool show_history = 0;
 bool make_dense   = 0;
 bool export_history = 0;
 bool create_history_file = 1;
 bool create_subdir  = 0;
 string in_positivityConstr;
 string in_GreenBasis;
 string label_eigen;
 string label_expect;
 string out_sdpa;
 string subdirSuffix = "";
 Vector<string> in_poly;
 OpSpec costOp_spec;
 Vector<OpSpec> eigen_spec;
 Vector<double> eigen_value;
 Vector<OpSpec> expect_spec;
 Vector<double> expect_value;
 int  prec = 12;
 bool prec_scientific = 0;
 bool prec_fixed = 0;
 size_t lineLen_max = 72;
 
 
 //=================================================================================
 //-[read in command line parameters]-----------------------------------------------
 if( argc<2) return 1;
 
 int argMode = 0;
 set<string> argModable;
 map<string,int> argID;
 argID["-in_positivityConstr"]	= 1;		//argModable.insert("-opLabels");
 argID["-in_GreenBasis"]	= 3;		//argModable.insert("-opLabels");
 argID["-set_in_poly"]		= 5;
 argID["-costOp"]		= 7;
 argID["-constrain_eigen"]	= 10;		argModable.insert("-constrain_eigen");
 argID["-label_eigen"]		= 11;
 argID["-constrain_expect"]	= 14;		argModable.insert("-constrain_expect");
 argID["-label_expect"]		= 15;
 argID["-show_history"]		= 20;
 argID["-make_dense"]		= 30;
 argID["-export_history"]	= 32;
 argID["-dont_create_history_file"] = 34;
 argID["-prec"]			= 40;
 argID["-prec_scientific"]	= 42;
 argID["-prec_fixed"]		= 44;
 argID["-out_sdpa"]		= 50;
 argID["-create_subdir"]	= 52;
 argID["-create_subdirSuffixed"]= 54;
 
 string curr_poly;
 for(int argn=1; argn<argc; argn++)
  {
  int aID;  string argStr( argv[argn] );
  if( argStr[0] == '-' && (argStr[1]<'0'||argStr[1]>'9') )
   { aID = argID[argStr];  argMode = 0;  if(argModable.count(argStr)) argMode = aID; }
  else { aID=argMode;  argn--; }
  switch( aID )
   {
   case   1:  in_positivityConstr = argv[argn+1];  argn++;  break;
   case   3:  in_GreenBasis       = argv[argn+1];  argn++;  break;
   case   5:  curr_poly           = argv[argn+1];  argn++;
              in_poly << curr_poly;  break;
   case   7:  costOp_spec  = OpSpec( curr_poly, argv[argn+1] );  argn++;  break;
   case  10:  eigen_spec  << OpSpec( curr_poly, argv[argn+1] );  argn++;
              eigen_value << boost::lexical_cast<double>( argv[argn+1] );  argn++;  break;
   case  11:  label_eigen = argv[argn+1];  argn++;  break;
   case  14:  expect_spec  << OpSpec( curr_poly, argv[argn+1] );  argn++;
              expect_value << boost::lexical_cast<double>( argv[argn+1] );  argn++;  break;
   case  15:  label_expect = argv[argn+1];  argn++;  break;
   case  20:  show_history   = 1;  break;
   case  30:  make_dense     = 1;  break;
   case  32:  export_history = 1;  break;
   case  34:  create_history_file = 0;  break;
   case  40:  prec = boost::lexical_cast<int>( argv[argn+1] );  argn++;  break;
   case  42:  prec_scientific = 1;  break;
   case  44:  prec_fixed      = 1;  break;
   case  50:  out_sdpa = argv[argn+1];  argn++;  break;
   case  52:  create_subdir   = 1;  break;
   case  54:  create_subdir   = 1;
              subdirSuffix = argv[argn+1];  argn++;  break;
   default :  PANIC("Wrong command line arguments!")(argStr);  break;
   }
  }
 
 PRECONDITION( fileExists(in_positivityConstr) )( in_positivityConstr );
 if( in_GreenBasis.size() )
  { PRECONDITION( fileExists(in_GreenBasis) )( in_GreenBasis ); }
 for(size_t n = 0; n < in_poly.size(); n++)
  { PRECONDITION( fileExists(in_poly[n]) )( in_poly[n] ); }
 PRECONDITION( !out_sdpa.size() || !create_subdir );
 
 
 cout << "\n===================== read PositivityConstraints, GreenBasis, and cost operator ==\n";
 // ---- PositivityConstraints ----
 History hist_posConstr;
 Vector<PositivityConstraintsF> constr;
 read( in_positivityConstr, hist_posConstr, constr );
 string file_GreenBasis = find_file( hist_posConstr[-1].get_param_s("file_GreenBasis"),
                                     in_positivityConstr );
 
 if( show_history )  { TRACE( hist_posConstr ); }

 // ---- GreenBasis ----
 History    hist_G;
 GreenBasis G_basis;
 LabelList<string> qnNames;
 load_cmp( in_GreenBasis, file_GreenBasis, hist_G, G_basis, qnNames );
 
 if( show_history )  { TRACE( hist_G ); }
 
 // ---- cost operator and constraint operators ----
 History hist_poly;
 PolynomialQn costOp;
 Vector<PolynomialQn> eigen_op(  eigen_spec.size() );
 Vector<PolynomialQn> expect_op( expect_spec.size() );
 {
 map<string,LabelledList<PolynomialQn> > list;
 for(size_t n = 0; n < in_poly.size(); n++)
  {
  LabelList<string> qnNames2;
  History hist_p;
  read( in_poly[n], hist_p, list[in_poly[n]], qnNames2 );
  if( in_poly[n] == costOp_spec.first )  hist_poly = hist_p;
  if( show_history )  { TRACE( in_poly[n] )( hist_p ); }
  CHECK_EQUAL( qnNames2, qnNames );
  CHECK_EQUAL( (int)hist_p[-1].get_info_d("algebraF"), G_basis.algebraF() );
  }
 PRECONDITION( list[costOp_spec.first].contains( costOp_spec.second ) );
 costOp = list[costOp_spec.first][ costOp_spec.second ];
 for(size_t n = 0; n < eigen_op.size(); n++)
  {
  PRECONDITION( list[eigen_spec[n].first].contains( eigen_spec[n].second ) );
  eigen_op[n] = list[eigen_spec[n].first][ eigen_spec[n].second ];
  }
 for(size_t n = 0; n < expect_op.size(); n++)
  {
  PRECONDITION( list[expect_spec[n].first].contains( expect_spec[n].second ) );
  expect_op[n] = list[expect_spec[n].first][ expect_spec[n].second ];
  }
 }
 
 
 cout << "\n======================================================================= History ==\n";
 History hist( hist_posConstr );
 hist.push_back( hist_poly[-1] );
 
 History::Action a;
 a.info_s["program"]      = program;
 a.info_s["command_line"] = commandLine( argc, argv );
 a.info_s["cwd"]          = current_workingPath();
 a.info_s["date"]         = dateTimeString();
 a.param_s["file_PositivityConstr"] = in_positivityConstr;
 a.param_s["file_GreenBasis"] = in_GreenBasis;
 a.param_s["costOp_file"]     = costOp_spec.first;
 a.param_s["costOp_label"]    = costOp_spec.second;
 a.param_d["create_subdir"] = create_subdir;
 a.param_d["make_dense"]    = make_dense;
 for(size_t n = 0; n < eigen_op.size(); n++)
  a.param_d["eigen_" +eigen_spec[n].second ] = eigen_value[n];
 for(size_t n = 0; n < expect_op.size(); n++)
  a.param_d["expect_"+expect_spec[n].second] = expect_value[n];
 hist.push_back(a);
 
 if( eigen_op.size() && !label_eigen.size() )
  {
  label_eigen = "eigen-";
  for(size_t n = 0; n < eigen_op.size(); n++)
   label_eigen += eigen_spec[n].second + boost::lexical_cast<string>( eigen_value[n] );
  }
 if( label_eigen.size() )
  hist[-1].info_s["label_PositivityConstraints_eigen"] = label_eigen;
 
 if( expect_op.size() && !label_expect.size() )
  {
  label_expect = "expect-";
  for(size_t n = 0; n < expect_op.size(); n++)
   label_expect += expect_spec[n].second + boost::lexical_cast<string>( expect_value[n] );
  }
 if( label_expect.size() )
 hist[-1].info_s["label_PositivityConstraints_expect"] = label_expect;
 
 
 cout << "\n========================================================== expand cost operator ==\n";
 costOp = normalOrdered( costOp, G_basis.algebraF() );
 CHECK_EQUAL( costOp, herm(costOp) );
 GreenFunctionDualHerm costV = greenFunctionDualHerm( G_basis.expand_dual( costOp ) );
 
 hist[-1].info_d["costOp_const"] = costV.first[0];
 
 
 if( eigen_op.size() )
  {
  cout << "\n========================================================= add eigen constriants ==\n";
  size_t NoBlocks = constr.size();
  resize( constr, NoBlocks + eigen_op.size(), true );
  for(size_t n = 0; n < eigen_op.size(); n++)
   {
   MonomialQn Id;
   PolynomialQn op = eigen_op[n] - eigen_value[n]*Id;
   op = normalOrdered( (-1)*op*op, G_basis.algebraF() );
   CHECK_EQUAL( op, herm(op) );
   constr[ NoBlocks + n ] = PositivityConstraintsF::from_observable( op, G_basis );
   }
  }
 
 
 if( expect_op.size() )
  {
  cout << "\n================================================ add expectation value constriants ==\n";
  size_t NoBlocks = constr.size();
  resize( constr, NoBlocks + 2*expect_op.size(), true );
  for(size_t n = 0; n < expect_op.size(); n++)
   {
   MonomialQn Id;
   PolynomialQn op = expect_op[n] - expect_value[n]*Id;
   op = normalOrdered( op, G_basis.algebraF() );
   CHECK_EQUAL( op, herm(op) );
   constr[ NoBlocks + 2*n ]   = PositivityConstraintsF::from_observable( op, G_basis );
   constr[ NoBlocks + 2*n+1 ] = PositivityConstraintsF::from_observable( (-1)*op, G_basis );
   }
  }
 
 
 cout << "\n==================================================== export to sdpa (text) file ==\n";
 if( !out_sdpa.size() )
  {
  out_sdpa = suffix_stripped( pathName_stripped( in_positivityConstr ), ".Constr.bin" );
  if( label_eigen.size()  )  out_sdpa += "_" + label_eigen;
  if( label_expect.size() )  out_sdpa += "_" + label_expect;
  out_sdpa += "_" + hist_poly[-1].get_info_s("param_Hamiltonian");
  string SubdirSuffix = ".export" + subdirSuffix;
  if( create_subdir )  system( (string("mkdir ")+out_sdpa+SubdirSuffix).c_str() );
  out_sdpa += string( create_subdir? SubdirSuffix+"/prob" : "" )
                  + ( make_dense?    ".sdpa.dat" : ".sdpa.dat-s" );
  }
 
 ofstream txt;
 txt.open( out_sdpa.c_str(), ios_base::ate );
 txt << setprecision(prec);
 if( prec_scientific )  txt << scientific;
 if( prec_fixed )       txt << fixed;
 
 string header = "command line: " + commandLine( argc, argv ) + "\n";
 if( export_history )
  {
  ostringstream oss;
  oss << hist;
  header += oss.str();
  }
 txt << "* " + boost::algorithm::replace_all_copy( wrap_lines( header, lineLen_max ), "\n", "\n* " );
 
 sdp_gs::Indices blockSize( constr.size() );
 for(size_t n = 0; n < constr.size(); n++)  blockSize[n] = constr[n].dim();
 
 
 txt << endl << G_basis.size()-1 << " = G_basis.size()-1";
 txt << endl << constr.size()    << " = constr.size()";
 txt << endl << blockSize        << " = (constr[n].dim())_n";
 
 txt << endl << direct_sum( costV.first[ range(1,costV.first.size()) ], costV.second );
 
 CHECK_EQUAL( G_basis.elementsReal()[0], MonomialQn() );
 
 // Id
 txt << endl;
 for(size_t m = 0; m < constr.size(); m++)
  txt << endl << print_matrix( make_dense, 0, m, -constr[m].matricesReal()[0] );
 // constraint matrices coupling to Re(G)
 for(size_t n = 1; n < G_basis.sizeReal(); n++)
  {
  txt << endl;
  for(size_t m = 0; m < constr.size(); m++)
   txt << endl << print_matrix( make_dense, n, m, constr[m].matricesReal()[n] );
  }
 // constraint matrices coupling to Im(G)
 for(size_t n = 0; n < G_basis.sizeImag(); n++)
  {
  txt << endl;
  for(size_t m = 0; m < constr.size(); m++)
   txt << endl << print_matrix( make_dense, G_basis.sizeReal()+n, m,
                                constr[m].matricesImag()[n] );
  }
 
 txt << endl;
 txt.close();
 
 cout << " minimization problem exported to SDPA file " << out_sdpa << "\n";
 
 
 if( create_history_file )
  {
  cout << "\n========================================================= creating history file ==\n";
  string out_hist = out_sdpa.substr( 0, out_sdpa.rfind(".") ) + ".history.bin";
  write( out_hist, hist );
  
  cout << " History exported to the file " << out_hist << "\n";
  }
 
 
 cout<<endl<<"END"<<endl;
 #if !defined(NDEBUG)
 cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 }
