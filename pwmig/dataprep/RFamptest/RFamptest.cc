#include <string>
#include <iostream>
#include <list>
#include <vector>
#include "perf.h"
#include "seispp.h"
#include "Metadata.h"
#include "dbpp.h"
#include "ensemble.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "RFamptest db [-use_wfdisc]" <<endl;
    exit(-1);
}
MetadataList CreateMdlForThis()
{
    MetadataList mdl;
    Metadata_typedef t;
    t.tag=string("sta");
    t.mdt=MDstring;
    mdl.push_back(t);
    t.tag=string("datatype");
    mdl.push_back(t);
    t.tag=string("dir");
    mdl.push_back(t);
    t.tag=string("dfile");
    mdl.push_back(t);
    t.mdt=MDint;
    t.tag=string("foff");
    mdl.push_back(t);
    t.tag=string("nsamp");
    mdl.push_back(t);
    t.mdt=MDreal;
    t.tag=string("time");
    mdl.push_back(t);
    t.tag=string("endtime");
    mdl.push_back(t);
    t.tag=string("samprate");
    mdl.push_back(t);
    return(mdl);
}

bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc<2) usage();
    string dbname(argv[1]);
    int i;
    bool use_wfdisc(false);
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-use_wfdisc")
            use_wfdisc=true;
        else
            usage();
    }
    try{
        DatascopeHandle dbh(dbname,true);
        list<string> sort_keys;
        if(use_wfdisc)
        {
            dbh.lookup("wfdisc");
            sort_keys.push_back("time");
            sort_keys.push_back("sta");
            sort_keys.push_back("chan");
        }
        else
        {
            dbh.lookup("wfprocess");
            dbh.natural_join("sclink");
            dbh.natural_join("evlink");
            sort_keys.push_back("pwfid");
        }
        dbh.sort(sort_keys);
        AttributeMap am;
        MetadataList mdl;
        mdl=CreateMdlForThis();
        Metadata_typedef t;
        if(use_wfdisc)
        {
            t.tag=string("chan");
            t.mdt=MDstring;
            mdl.push_back(t);
        }
        else
        {
            t.tag=string("timetype");
            t.mdt=MDstring;
            mdl.push_back(t);
            t.mdt=MDint;
            t.tag=string("evid");
            mdl.push_back(t);
        }
        int nrows=dbh.number_tuples();
        for(i=0,dbh.rewind();i<nrows;++i,++dbh)
        {
            double rms;
            Metadata *mptr;
            Metadata md;
            /* This is an ugly consequence of multiple inheritance */
            if(use_wfdisc)
            {
                TimeSeries *d;
                d=new TimeSeries(dbh,mdl,am);
                rms=dnrm2(d->ns,&(d->s[0]),1);
                rms /= static_cast<double>(d->ns);
                mptr=dynamic_cast<Metadata *>(d);
                md=Metadata(*mptr);
                delete d;
            }
            else
            {
                ThreeComponentSeismogram *d;
                d=new ThreeComponentSeismogram(dbh,mdl,am);
                int nts=3*d->ns;
                rms=dnrm2(nts,d->u.get_address(0,0),1);
                rms /= static_cast<double>(nts);
                mptr=dynamic_cast<Metadata *>(d);
                md=Metadata(*mptr);
                delete d;
            }
            cout << md.get_string("sta")<<" ";
            if(use_wfdisc)
                cout << md.get_string("chan")<<" ";
            else
                cout << md.get_int("evid")<<" ";

            cout <<strtime(md.get_double("time"))<< " " << rms <<endl;
        }
    }catch(SeisppError& serr)
    {
        serr.log_error();
        exit(-2);
    }
    catch(std::exception& stex)
    {
        cerr << stex.what()<<endl;
        exit(-3);
    }
}
