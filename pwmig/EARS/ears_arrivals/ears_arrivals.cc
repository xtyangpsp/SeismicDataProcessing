#include <stdlib.h>
#include <string>
#include "stock.h"
#include "seispp.h"
#include "dbpp.h"

using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "ears_arrivals dbin dbout [-lag lag]"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc<3) usage();
    string dbin(argv[1]);
    string dbout(argv[2]);
    double lag0(10.0);
    int i;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-lag") 
        {
            ++i;
            if(i>=argc) usage();
            sarg=string(argv[i]);
            lag0=atof(sarg.c_str());
        }
        else
            usage();
    }
    try {
        DatascopeHandle dbhi(dbin,true);
        dbhi.lookup(string("wfdisc"));
        string subss("chan=~/ITR.*/");
        dbhi.subset(subss);
        cout << "ears_arrivals subset view of wfdisc number of rows="<<dbhi.number_tuples()<<endl;
        DatascopeHandle dbho(dbout,false);
        dbho.lookup(string("arrival"));
        long nrows;
        nrows=dbhi.number_tuples();
        long row;
        for(dbhi.rewind(),row=0;row<nrows;++row,++dbhi)
        {
            string sta=dbhi.get_string("sta");
            string chan=dbhi.get_string("chan");
            // Very specialized.  Changes the first 3 chars to BHZ 
            chan[0]='B';
            chan[1]='H';
            chan[2]='Z';
            double time0=dbhi.get_double("time");
            double tp=time0+lag0;
            dbho.append();
            dbho.put("sta",sta);
            dbho.put("chan",chan);
            dbho.put("time",tp);
            dbho.put("jdate",yearday(tp));
            string iphase("P");
            dbho.put("iphase",iphase);
            double deltim(0.1);
            dbho.put("deltim",deltim);
        }
    } catch (SeisppError& serr)
    {
        serr.log_error();
        exit(-1);
    }
}



