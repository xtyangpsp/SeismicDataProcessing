#include <stdio.h>
#include <string>
#include <iostream>
#include "db.h"
using namespace std;
void usage()
{
    cerr << "earskeepers dbin dbout"<<endl
        << "dbin wfdisc with sac files having user8 set will be copied to dbout"
        <<endl;
    exit(-1);
}
int main(int argc, char **argv)
{
    ios::sync_with_stdio();
    if(argc!=3) usage();
    long nrec;
    Dbptr db;
    int iret;
    iret=dbopen(argv[1],"r",&db);
    db=dblookup(db,NULL,"wfdisc",NULL,NULL);
    if(iret==dbINVALID)
    {
        cerr << "Open failed on input db "<< argv[1]<<endl;
        usage();
    }
    Dbptr dbout;
    iret=dbopen(argv[2],"r+",&dbout);
    if(iret==dbINVALID)
    {
        cerr << "Open failed on output db "<< argv[2]<<endl;
        usage();
    }
    cout << "Reading from wfdisc of database "<<argv[1]<<endl
        << "Writing results to database "<<argv[2]<<endl;
    dbout=dblookup(dbout,NULL,"wfdisc",NULL,NULL);
    dbquery(db,dbRECORD_COUNT,static_cast<void *>(&nrec));
    long nkeepers(0);
    for(db.record=0;db.record<nrec;++db.record)
    {
        char dir[128],dfile[40];
        char fullrecord[512];
        iret=dbgetv(db,0,"dir",dir,"dfile",dfile,NULL);
        if(iret==dbINVALID)
        {
            cerr << "dbgetv error reading wfdisc record "
                << db.record<<endl
                << "Fix problem and try again"<<endl;
            exit(-1);
        }
        string fname=string(dir)+"/"+string(dfile);
        /* We could do this fancy but use a more fragile
           and less general approach here.  seek and extract
           only the one key header field (user8).  If set
           true copy the wfdisc records.  if false (0) skip. */
        FILE *fp;
        fp=fopen(fname.c_str(),"r");
        long user8_offset(192);
        fseek(fp,user8_offset,SEEK_SET);
        float user8;
        size_t count;
        count=fread(&user8,sizeof(float),1,fp);
        if(user8>0.0)
        {
            cout << "Copying record for sac file="<<fname<<endl;
            dbget(db,fullrecord);
            dbadd(dbout,fullrecord);
            ++nkeepers;
        }
    }
    cout << "Finished:   Wrote "<<nkeepers<<" records to output wfdisc"<<endl
        << "Input db number of rows = "<<nrec<<endl;
}

