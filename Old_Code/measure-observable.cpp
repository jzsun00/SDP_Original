#include "io/operator00_io.h"
#include "io/symmetries00_io.h"
#include "io/GreenFunction00_io.h"

#include "fileioL00.h"
#include "linalg/MatrixTableL00_io.h"
#include "linalg/TensorL02.h"
#include "parseMathL03.h"
#include "labelledValueL02.h"

#include "boost/algorithm/string/trim.hpp"
#include "boost/algorithm/string/split.hpp"


using namespace std;
using namespace LinearAlgebra;



//##################################################################################################

// TODO: For strange reasons we need to make sure explicitely that boost does not play around with Vector<ParseMath::Formula> when including boost/algorithm/string/trim.hpp.
namespace boost
{
template<>
struct range_difference<Vector<ParseMath::Formula> >
 {
 };
}


//--------------------------------------------------------------------------[ is_switch ]-
bool is_switch(string const& argStr)
 {
 return( argStr[0] == '-' && (argStr[1]<'0'||argStr[1]>'9') );
 }


//---------------------------------------------------------------------[ read_strVector ]-
Vector<string> read_strVector(string const& txt)
 {
 Vector<string> ret;
 size_t pos  = txt.find("(");
 size_t rpos = txt.rfind(")");
 if( pos == string::npos || rpos == string::npos )  return ret;
 string Txt = txt.substr( pos+1, rpos-(pos+1) ) + ",";
 pos = 0;
 do
  {
  size_t pos0 = pos;
  pos = Txt.find( ",", pos0 ) + 1;
  ret << boost::trim_copy( Txt.substr( pos0, pos-pos0-1 ) );
  }
 while( pos < Txt.size() );
 return ret;
 }


//----------------------------------------------------------------------[ OperatorDescr ]-
struct OperatorDescr
 {
 OperatorDescr()  {}
 OperatorDescr(string const& i_label_in, string const& i_label)
  : label_in(i_label_in), label(i_label)  {}
 OperatorDescr(string const& i_label_in, string const& i_label, Vector<string> const& i_shift)
  : label_in(i_label_in), label(i_label), shift(i_shift.size())
  { for(size_t n = 0; n < shift.size(); n++)  shift[n] = ParseMath::parse( i_shift[n] ); }
 
 string label_in;
 string label;
 Vector<ParseMath::Formula> shift;
 };


//---------------------------------------------------------------------[ Measure_opProd ]-
struct Measure_opProd
 {
 Measure_opProd()  {}
 Measure_opProd(string const& i_label, string const& opProd_str) : label(i_label)
  {
  vector<string> splitVec;
  boost::split( splitVec, opProd_str, boost::is_any_of("-*") );
  opLabel = splitVec;
  }
 
 string label;
 Vector<string> opLabel;
 };



//__________________________________________________________________________________________________
//########################################################################################### MAIN #
int main(int argc, char* argv[])
 {
 bool silent=0;
 bool supersilent=0;
 for(int argn=1; argn<argc; argn++)
  {
  if( string(argv[argn]) == "-silent" )       silent = 1;
  if( string(argv[argn]) == "-supersilent" )  silent = supersilent = 1;
  }
 
 #if !defined(NDEBUG)
 if(!silent)  cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 
 string program = "measure-observable";
 
 if(!silent)
  {
  cout<<" >>> Measure observable(s) for a given Green's function <<< 29.08.2011\n";
  cout<<" usage:  "<< program << endl;
  cout<<"		 -in_G <*.G.bin file>  (-in_GreenBasis <*.Gbasis.bin file>)  (-ver <requested program version =1>)\n";
  cout<<"		(-set_in_poly <*.PolyQn.list.bin file>)  (-show_history)  (-show_params)  (-show_opLabels)\n";
  cout<<"		(-op <operator label from the current poly file> ...)  (-op_renamed <operator label from the current poly file> <new label> ...)\n";
  cout<<"		(-op_shifted <operator label> <latt. shift vector> <new label> ...)\n";
  cout<<"		(-add_loop <variabel label> <N> ...)  (-add_var <variabel label> <formula> ...)\n";
  cout<<"		(-precondition <formula>)\n";
  cout<<"		(-measure_prod <result label> <op labels seperated by '*'> ...)\n";
  cout<<"		(-col <col label> <formula> ...)\n";
  cout<<"		(-out_table <*.Table.bin .. in case you wanna store the result>)\n";
  cout<<"		(-dont_disp_table)  (-disp_listRow)  (-silent)  (-supersilent)\n";
  cout<<"		(-sep <column separator ='\\t'>)  (-prec <number of digits =12>)  (-prec_scientific)  (-prec_fixed)\n";
  cout<<" You can call '-set_in_poly' several times. Each time, use '-op' or '-op_rename' to specify which operators\n";
  cout<<"  from the current file are to be used.\n";
  cout<<" Loop variables introduced via '-add_loop' run from 0 to N-1.\n";
  cout<<" Variable evaluation sequence: (a) parameters from the histories (A: *.G.bin, B: *.Gbasis.bin, C,..: *.PolyQn.list.bin),\n";
  cout<<"  (b) loop variables, (c) add. variables, (d) precondition, (e) measurement results.\n";
  cout<<endl;
  }
 cout << setprecision(12);
 
 bool show_history = 0;
 bool show_params  = 0;
 bool show_opLabels = 0;
 bool disp_table   = 1;
 bool disp_listRow = 0;
 string sep = "\t";
 string in_G;
 string in_GreenBasis;
 Vector<string>          in_poly;
 Vector<Vector<OperatorDescr> > in_op;
 Vector<string>          loop_label;
 Vector<size_t>          loop_dim;
 Vector<string>              addVar_label;
 Vector<ParseMath::Formula>  addVar;
 ParseMath::Formula          precondition = 1;
 Vector<Measure_opProd>      measure_prod;
 Vector<string>              colLabel;
 Vector<ParseMath::Formula>  col;
 string out_table;
 
 
 //=================================================================================
 //-[read in command line parameters]-----------------------------------------------
 if( argc<2) return 1;
 
 int argMode = 0;
 set<string> argModable;
 map<string,int> argID;
 argID["-ver"]			= 1;
 argID["-in_G"]			= 3;
 argID["-in_GreenBasis"]	= 5;
 argID["-show_history"]		= 10;
 argID["-show_params"]		= 12;
 argID["-show_opLabels"]	= 14;
 argID["-set_in_poly"]		= 20;
 argID["-op"]			= 22;		argModable.insert("-op");
 argID["-op_renamed"]		= 24;		argModable.insert("-op_renamed");
 argID["-op_shifted"]		= 26;		argModable.insert("-op_shifted");
 argID["-add_loop"]		= 30;		argModable.insert("-add_loop");
 argID["-add_var"]		= 40;		argModable.insert("-add_var");
 argID["-precondition"]		= 46;
 argID["-measure_prod"]		= 60;		argModable.insert("-measure_prod");
 argID["-col"]			= 70;		argModable.insert("-col");
 argID["-out_table"]		= 80;
 argID["-dont_disp_table"]	= 135;
 argID["-disp_listRow"]		= 136;
 argID["-silent"]		= 137;
 argID["-supersilent"]		= 138;
 argID["-sep"]			= 139;
 argID["-prec"]			= 140;
 argID["-prec_scientific"]	= 142;
 argID["-prec_fixed"]		= 144;
 
 for(int argn=1; argn<argc; argn++)
  {
  int aID;  string argStr( argv[argn] );
  if( is_switch(argStr) )
   { aID = argID[argStr];  argMode = 0;  if(argModable.count(argStr)) argMode = aID; }
  else { aID=argMode;  argn--; }
  switch( aID )
   {
   case   1:  CHECK_EQUAL( boost::lexical_cast<int>( argv[argn+1] ), 1 );  argn++;  break;
   case   3:  in_G          = argv[argn+1];  argn++;  break;
   case   5:  in_GreenBasis = argv[argn+1];  argn++;  break;
   case  10:  show_history  = 1;  break;
   case  12:  show_params   = 1;  break;
   case  14:  show_opLabels = 1;  break;
   case  20:  in_poly << argv[argn+1];  argn++;
              in_op   << Vector<OperatorDescr>();  break;
   case  22:  in_op[ in_op.size()-1 ] << OperatorDescr( argv[argn+1], argv[argn+1] );  argn++;   break;
   case  24:  in_op[ in_op.size()-1 ] << OperatorDescr( argv[argn+1], argv[argn+2] );  argn+=2;  break;
   case  26:  in_op[ in_op.size()-1 ] << OperatorDescr( argv[argn+1], argv[argn+3],
                                                        read_strVector(argv[argn+2]) );  argn+=3;  break;
   case  30:  loop_label << argv[argn+1];
              loop_dim   << boost::lexical_cast<size_t>( argv[argn+2] );  argn+=2;  break;
   case  40:  addVar_label << argv[argn+1];
              addVar       << ParseMath::parse( argv[argn+2] );  argn+=2;  break;
   case  46:  precondition = ParseMath::parse( argv[argn+1] );  argn++;  break;
   case  60:  measure_prod << Measure_opProd( argv[argn+1], argv[argn+2] );  argn+=2;  break;
   case  70:  colLabel << argv[argn+1];  argn++;
              col << ParseMath::parse( argv[argn+1] );  argn++;  break;
   case  80:  out_table = argv[argn+1];  argn++;  break;
   
   case 135:  disp_table   = 0;  break;
   case 136:  disp_listRow = 1;  break;
   case 137:  PRECONDITION( silent );  break;
   case 138:  PRECONDITION( silent && supersilent );  break;
   case 139:  sep = argv[argn+1];  argn++;  break;
   case 140:  cout << setprecision( boost::lexical_cast<int>( argv[argn+1] ) );  argn++;  break;
   case 142:  cout << scientific;  break;
   case 144:  cout << fixed;       break;
   default :  PANIC("Wrong command line arguments!")(argStr);  break;
   }
  }
 
 PRECONDITION( fileExists(in_G) )( in_G );
 if( in_GreenBasis.size() )
  { PRECONDITION( fileExists(in_GreenBasis) )( in_GreenBasis ); }
 for(size_t n = 0; n < in_poly.size(); n++)
  { PRECONDITION( fileExists( in_poly[n] ) )( in_poly[n] ); }
 
 
 if(!silent)  cout << "\n================================= read GreenFunction, GreenBasis, and operators ==\n";
 map<string,double> param;
 map<string,string> paramString;	// not used
 
 // ---- GreenFunction ---- 
 History hist_G;
 GreenFunction G;
 read( in_G, hist_G, G );
 string file_GreenBasis = find_file( hist_G[-1].get_param_s("file_GreenBasis"), in_G );
 hist_G.export_all( paramString, param, "A" );
 
 if( show_history )  { TRACE( hist_G ); }
 
 // ---- GreenBasis ----
 History    hist_Gb;
 GreenBasis G_basis;
 LabelList<string> qnNames;
 load_cmp( in_GreenBasis, file_GreenBasis, hist_Gb, G_basis, qnNames );
 hist_Gb.export_all( paramString, param, "B" );
 
 if( show_history )  { TRACE( hist_Gb ); }
 
 Lattice const& latt( G_basis.symmetries().lattice().symm );
 sdp_gs::Indices const& coord_idcs( G_basis.symmetries().lattice().idcs );
 
 // ----  Polynomials ---- 
 typedef LabelledList<Vector<ParseMath::Formula> > LabelledShiftList;
 History hist_poly;
 LabelledList<PolynomialQn> ops;
 LabelledShiftList          op_shifts;
 for(size_t n = 0; n < in_poly.size(); n++)
  {
  LabelledList<PolynomialQn> list;
  LabelList<string> qnNames2;
  read( in_poly[n], hist_poly, list, qnNames2 );
  if( show_history )   { TRACE( in_poly[n] )( hist_poly ); }
  if( show_opLabels )  { TRACE( in_poly[n] )( list.labels() ); }
  CHECK_EQUAL( qnNames2, qnNames )( in_poly[n] );
  CHECK_EQUAL( (int)hist_poly[-1].get_info_d("algebraF"), G_basis.algebraF() );
  hist_poly.export_all( paramString, param, boost::lexical_cast<string>(char('C'+n)) );
  
  for(size_t nOp = 0; nOp < in_op[n].size(); nOp++)
   {
   PRECONDITION( list.contains( in_op[n][nOp].label_in ) )( in_op[n][nOp].label_in );
   ops[       in_op[n][nOp].label ] = list.find( in_op[n][nOp].label_in )->second;
   op_shifts[ in_op[n][nOp].label ] = in_op[n][nOp].shift;
   }
  }
 
 // Parameters
 if( show_params )  { TRACE( param );  return -1; }
 ParseMath::VariableMap varMap( param );
 
 
 if(!silent)  cout << "\n======================================================================= History ==\n";
 History hist( hist_G );
 
 History::Action a;
 a.info_s["program"]      = program;
 a.info_s["command_line"] = commandLine( argc, argv );
 a.info_s["cwd"]          = current_workingPath();
 a.info_s["date"]         = dateTimeString();
 a.param_s["file_GreenFunction"] = in_G;
 hist.push_back(a);
 
 
 if(!silent)  cout << "\n========================================================= pre-evaluate formulas ==\n";
 precondition = ParseMath::eval( precondition, varMap );
 for(size_t n = 0; n < addVar.size(); n++)  addVar[n] = ParseMath::eval( addVar[n], varMap );
 for(size_t n = 0; n < col.size();    n++)  col[n]    = ParseMath::eval( col[n],    varMap );
 for(LabelledShiftList::iterator i = op_shifts.begin(); i != op_shifts.end(); ++i)
  for(size_t m = 0; m < i->second.size(); m++)
   i->second[m] = ParseMath::eval( i->second[m], varMap );
 
 
 if(!silent)  cout << "\n============================================== do loop and evaluate observables ==\n";
 set<Vector<double> > dat;
 
 TensorIndex ti( loop_dim );
 for(TensorIndexState i = ti.begin(); i.is_valid(); ++i)
  {
  ParseMath::VariableMap val;
  
  // evaluate add. variables
  for(size_t n = 0; n < ti.rank(); n++)
   val[ loop_label[n] ] = i[n];
  for(size_t n = 0; n < addVar.size(); n++)
   val[ addVar_label[n] ] = ParseMath::eval( addVar[n], val );
  
  if( 0 != ParseMath::eval_double( precondition, val ) )
   {
   // get (and shift) measurement operators
   LabelledList<PolynomialQn> opShifted;
   for(LabelledList<PolynomialQn>::const_iterator i = ops.begin(); i != ops.end(); ++i)
    {
    PolynomialQn Op;
    if( op_shifts[i->first].size() )
     {
     sdp_gs::PointI shift( op_shifts[i->first].size() );
     for(size_t d = 0; d < shift.size(); d++)
      {
      double dx = ParseMath::eval_double( op_shifts[i->first][d], val );
      CHECK_EQUAL( dx, (int)dx );
      shift[d] = (int)dx;
      }
     Op = translated( i->second, coord_idcs, latt, shift );
     }
    else  Op = i->second;
    opShifted[i->first] = Op;
    }
   
   // evaluate measurement operator products
   for(size_t n = 0; n < measure_prod.size(); n++)
    {
    if( !measure_prod[n].opLabel.size() )
     val[ measure_prod[n].label ] = 0;
    else
     {
     PolynomialQn Op = opShifted[ measure_prod[n].opLabel[0] ];
     for(size_t nOp = 1; nOp < measure_prod[n].opLabel.size(); nOp++)
      Op = Op*opShifted[ measure_prod[n].opLabel[nOp] ];
     val[ measure_prod[n].label ] = G_basis.eval( Op, G );
     }
    }
  
   // evaluate table row
   Vector<double> row( col.size(), 3e8 );
   for(size_t c = 0; c < col.size(); c++)
    {
    row[c] = ParseMath::eval_double( col[c], val );
    if( val.find( colLabel[c] ) == val.end() )
     val[ colLabel[c] ] = row[c];
    }
   dat.insert( row );
   }
  }
 
 
 if(!silent)  cout << "\n============================== build output MatrixTable and display if required ==\n";
 MatrixTable::labelList_type label2;
 for(size_t c = 0; c < col.size(); c++)
  label2.push_back( colLabel[c] );
 
 MatrixTable table( dat.size(), label2 );
 
 if( disp_table && !supersilent )
  {
  cout << "# " << commandLine( argc, argv );
  cout << "\n#\n# ";
  for(size_t c = 0; c < col.size(); c++)  cout << colLabel[c] << ((c+1)<col.size()? sep:"\n");
  }
 
 size_t row = 0;
 for(set<Vector<double> >::const_iterator it = dat.begin(); it != dat.end(); ++it, ++row)
  {
  if( disp_table )
   {
   if( disp_listRow )  cout << row << (col.size()? sep:"\n");
   for(size_t c = 0; c < col.size(); c++)
    cout << (*it)[c] << ((c+1)<col.size()? sep:"\n");
   }
  table.row( row ) = *it;
  }
 
 
 if( out_table.size() )
  {
  if(!silent)  cout << "\n============================================================= store MatrixTable ==\n";
  store( out_table, hist, table );
  }
 
 
 if(!silent)  cout<<endl<<"END"<<endl;
 #if !defined(NDEBUG)
 if(!silent)  cout<<"\n DON'T FORGET TO COMPILE WITH PRODUCTION-FLAGS, WHEN YOU'RE NOT DEBUGGING!\n\n";
 #endif
 }
