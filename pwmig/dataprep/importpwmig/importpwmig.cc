#include <string>
#include <list>
#include "gclgrid.h"
#include "seispp.h"
#include "dbpp.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "importpwmig [-vector|-scalar -a -db dbname] file1 file2 ... filen"<<endl
        <<  "  only pwmig output fields can be imported"<<endl
        <<  "  Use -vector migration output (default)and -scalar for coherence data"<<endl
        <<  "  dbname default is imports"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc,char **argv)
{
    if(argc<3) usage();
    bool isvector(true);
    string dbname("imports");
    string outdir("pwmigdata");
    string nulldir("");
    bool append_mode(false);
    list<string> fieldlist;
    /* First parse the arg list to build the list of field names to import*/
    int i;
    for(i=1;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-db")
        {
            ++i;
            if(i>=argc) usage();
            dbname=string(argv[i]);
        }
        else if(sarg=="-vector")
            isvector=true;
        else if(sarg=="-scalar")
            isvector=false;
        else if(sarg=="-a")
            append_mode=false;
        else
            fieldlist.push_back(sarg);
    }
    try {
        DatascopeHandle dbh(dbname,false);
        if(isvector)
        {
		    list<string>::iterator fptr;
		    for(fptr=fieldlist.begin(),i=0;fptr!=fieldlist.end();++fptr,++i)
		    {
		        try {
		            GCLvectorfield3d f(*fptr);
                            if( (i>0) || append_mode)
                                f.save(dbh,nulldir,outdir,*fptr,*fptr);
                            else
		                f.save(dbh,outdir,outdir,*fptr,*fptr);
		        }catch(GCLgridError& gerr)
		        {
		            cerr << "Problems with input file "<<*fptr<<endl
		                << "Error message:  "<<gerr.what()<<endl
		                << "Data from that file not imported"<<endl;
		        }
		    }
        }
        else
        {
		    list<string>::iterator fptr;
		    for(fptr=fieldlist.begin(),i=0;fptr!=fieldlist.end();++fptr,++i)
		    {
		        try {
		            GCLscalarfield3d f(*fptr);
                            if( (i>0) || append_mode)
                                f.save(dbh,nulldir,outdir,*fptr,*fptr);
                            else
		                f.save(dbh,outdir,outdir,*fptr,*fptr);
		        }catch(GCLgridError& gerr)
		        {
		            cerr << "Problems with input file "<<*fptr<<endl
		                << "Error message:  "<<gerr.what()<<endl
		                << "Data from that file not imported"<<endl;
		        }
		    }
        }
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(GCLgridError& gerr)
    {
        cerr << "Fatal Error:  "<<gerr.what()<<endl;
    }
}
