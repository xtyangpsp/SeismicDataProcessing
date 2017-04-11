#include <stdio.h>
#include "coords.h"
#include "seispp.h"
#include "EventCatalog.h"
#include "dbpp.h"
#include "SEEDStaGeom.h"
using namespace std;
using namespace SEISPP;
class SEEDSiteTable
{
public:
    SEEDSiteTable();
    /* Add if not yet define.  Should return 0 if added,
       1 if already defined, and -1 if already defined with
       different coordinates.  -1 case needs to be trapped for
       interactive error handling. */
    int add(SEEDStaGeom& sta);
    /* Unlike add this will force replacing if a duplicate. 
       Will also add if not yet defined */
    void replace(SEEDStaGeom& sta);
    /* Find station with net==n and sta==s */
    SEEDStaGeom find(string n, string s);
    /* Returns count of number of stations with name s*/
    int count(string s);
    /* Returns a list of net names that have station name s*/
    list<string> nets(string s);
    /* When an add fails we cache the duplicate found in
       the existing table.   Return the duplicate found.
       This is not at all a generic interface, but do not
       consider this object to be portable but a custom
       interface for this program only. */
    SEEDStaGeom last_duplicate(){
        return(duplicate);
    };
    void operator ++();
    void rewind();
    int size();
    SEEDStaGeom current();
private:
    map<string,SEEDStaGeom> sitetable;
    map<string,SEEDStaGeom>::iterator current_row;
    /* This is a cross reference with sta as key and net as
       value.  multimap count method is used to detect same name
       in different networks */
    multimap<string,string> xref;
    SEEDStaGeom duplicate;
    string key(string net, string sta);
};
SEEDSiteTable::SEEDSiteTable()
{
}
int SEEDSiteTable::add(SEEDStaGeom& s)
{
    int retcode;
    string testkey=this->key(s.net,s.sta);
    map<string,SEEDStaGeom>::iterator it,ite;
    ite=sitetable.end();
    it=sitetable.find(testkey);
    if(it==ite)
    {
        sitetable.insert(pair<string,SEEDStaGeom>(testkey,s));
        xref.insert(pair<string,string>(s.sta,s.net));
        retcode=0;
    }
    else
    {
        duplicate=(*it).second;
        if(((*it).second)==s)
            retcode=1;
        else
        {
            retcode=-1;
        }
    }
    return(retcode);
}
void SEEDSiteTable::replace(SEEDStaGeom& sta)
{
    /* I am not sure this will work if sta was not previously 
       defined.  For this program this won't happen, but don't 
       cut and paste this code without verifying this would
       create a seg fault bug */
    string thiskey=this->key(sta.net,sta.sta);
    map<string,SEEDStaGeom>::iterator it,ite;
    ite=sitetable.end();
    it=sitetable.find(thiskey);
    if(it==ite)
    {
        sitetable.insert(pair<string,SEEDStaGeom>(thiskey,sta));
        xref.insert(pair<string,string>(sta.sta,sta.net));
    }
    else
    {
        (*it).second=sta;
    }
    /* No need to touch xref multimap as this should only be 
       called if net,sta is already set for this station */
}
int SEEDSiteTable::count(string s)
{
    size_t n=xref.count(s);
    return(n);
}
list<string> SEEDSiteTable::nets(string s)
{
    list<string> result;
    multimap<string,string>::iterator it;
    pair<multimap<string,string>::iterator,multimap<string,string>::iterator> ret;
    ret=xref.equal_range(s);
    for(it=ret.first;it!=ret.second;++it)
        result.push_back((*it).second);
    return(result);
}
void SEEDSiteTable::operator++()
{
    if(current_row==sitetable.end())
        throw SeisppError(string("SEEDSiteTable ++operator:  ")
               +" tried to iterate past end of table");
    ++current_row;
}
void SEEDSiteTable::rewind()
{
    current_row=sitetable.begin();
}
SEEDStaGeom SEEDSiteTable::current()
{
    return((*current_row).second);
}
int SEEDSiteTable::size()
{
    return(sitetable.size());
}
string  SEEDSiteTable::key(string n, string s)
{
    return(n+"_"+s);
}
/* VERY fragile method for getting netname.  This parses the 
   input database full path name (path) and extracts the net
   name from the lowest level directory name.  This depends totally
   on a naming convention derived from EARS original layout.  
   Just BE WARNED this requires coordination with shell scripts
   that assemble data directory tree.

Assumed form of path string is similar to this:
/Volumes/SeismicData2/EARS2012/event_gathers/wf/evid1/II/IIdb

where in this case the net name is II

Changed June 26, 2012 when the file structure changed.
Now this is passed the dir field from a wfdisc.  Split is
simple as it is drive by names like II.PFO.  Just split to
II and PFO
 */
string extract_net_string(string path)
{
    size_t left,right;
    /* old version of this code
    right=path.rfind("/");
    left=path.rfind("/",right-1);
    string net;
    // this assign assumes left points at the first / and
    // right is the position of the right hand / 
    net.assign(path,left+1,right-left-1);
    */
    /* early 2012 version was this 
    right=path.find(".");
    string net;
    net.assign(path,0,right);
    */
    // This is version for station gathers EARS2013
    right=path.rfind("/");
    left=path.rfind(".");
    string net;
    net.assign(path,right+1,left-right-1);
    return(net);
}
string extract_directory_name(string path)
{
    size_t right;
    right=path.rfind("/");
    string dir;
    dir.assign(path,0,right);
    return(dir);
}
int save_origin(DatascopeHandle& dbh, EventCatalog& ec)
{
    try {
        DatascopeHandle dbev(dbh);
        dbh.lookup("origin");
        dbev.lookup("event");
        /* These are intentionally not fatal because we can likely
           repair this kind of problem after the fact */
        if(dbh.number_tuples()>0) 
            cerr << "Warning origin table in output db is not empty"
                <<endl;
        if(dbev.number_tuples()>0) 
            cerr << "Warning event table in output db is not empty"
                <<endl;
        ec.rewind();
        int nhypos=ec.size();
        cout << "Starting to write event and origin tables"<<endl
            << "Number of unique hypocenters found="<<nhypos<<endl;
        Hypocenter h;
        for(int i=0;i<nhypos;++i,++ec)
        {
            h=ec.current();
            /* Use counter + 1 as evid.  nass and ndef are 
            assigned constnt, bogus values because they are
            now set as primary keys  */
            long evid=i+1;
            long nass(10),ndef(10);
            long rec=dbaddv(dbh.db,0,"lat",deg(h.lat),
                    "lon",deg(h.lon),
                    "depth",h.z,
                    "time",h.time,
                    "evid",evid,
                    "orid",evid,
                    "nass",nass,
                    "ndef",ndef,NULL);
            if(rec<0) cerr << "Warning:  dbaddv failed "
                << "saving origin table row for evid "
                << evid<<endl;
            rec=dbaddv(dbev.db,0,"evid",evid,
                    "prefor",evid,NULL);
            if(rec<0) cerr << "Warning:  dbaddv failed "
                << "saving event table row for evid "
                << evid<<endl;
        }
    }catch (...){throw;}
}
string seedname(string net, string sta)
{
    string longname;
    if(sta.size()<4)
        longname=net+"_"+sta;
    else
        longname=net+sta;
    return(longname);
}
/* This saves entries in site, sitechan, and snetsta.
 sitechan writes generic Z,N,E channel codes. Returns multimap of 
 stations with duplicate names and the names assigned them for
 final site and sitechan table.*/
multimap<string,string> save_site(DatascopeHandle dbh, SEEDSiteTable& t)
{
    multimap<string,string> dups;
    const string base_error("make_earsdb save_site procedure:  ");
    try {
        DatascopeHandle dbsite(dbh);
        DatascopeHandle dbschan(dbh);
        dbh.lookup("snetsta");
        dbsite.lookup("site");
        dbschan.lookup("sitechan");
        t.rewind();
        int nrow=t.size();
        cout << "Attempting to write data for "
            << nrow << "stations"<<endl;
        long chanid(1);
        for(int i=0;i<nrow;++i,++t,++chanid)
        {
            long rec;
            string sta;
            long ondate(1990001);
            SEEDStaGeom s=t.current();
            int n=t.count(s.sta);
            if(n>1)
            {
                sta=seedname(s.net,s.sta);
                dups.insert( pair<string,string>(s.sta,sta) );
            }
            else
                sta=s.sta;
            rec=dbaddv(dbsite.db,0,"sta",sta.c_str(),
                    "ondate",ondate,
                    "lat",s.lat,
                    "lon",s.lon,
                    "elev",s.elev,
                    NULL);
            if(rec<0) throw SeisppError(base_error
                    + "dbaddv failed writing site table");
            rec=dbaddv(dbh.db,0,"snet",s.net.c_str(),
                    "fsta",s.sta.c_str(),
                    "sta",sta.c_str(),NULL);
            if(rec<0) throw SeisppError(base_error
                    + "dbaddv failed write snetsta table");
            rec=dbaddv(dbschan.db,0,"sta",sta.c_str(),
                    "chan","E",
                    "chanid",chanid,
                    "ondate",ondate,
                    "ctype","c",
                    "hang",90.0,
                    "vang",90.0,
                    NULL);
            /* Test only first dbaddv to sitechan because can't 
               conceive how the 2nd and 3rd would fail if the first
               succeeds. */
            if(rec<0) throw SeisppError(base_error
                    + "dbaddv failed write sitechan table");
            ++chanid;
            rec=dbaddv(dbschan.db,0,"sta",sta.c_str(),
                    "chan","N",
                    "chanid",chanid,
                    "ondate",ondate,
                    "ctype","c",
                    "hang",0.0,
                    "vang",90.0,
                    NULL);
            ++chanid;
            rec=dbaddv(dbschan.db,0,"sta",sta.c_str(),
                    "chan","Z",
                    "chanid",chanid,
                    "ondate",ondate,
                    "ctype","c",
                    "hang",0.0,
                    "vang",0.0,
                    NULL);
        }
        return(dups);
    }catch (...){throw;}
}
bool coords_bad(SEEDStaGeom& sta)
{
    if(sta.lat<-90.0 || sta.lat>90.0) return(true);
    if(sta.lon<-180 || sta.lon>180.0) return(true);
    return(false);
}
void usage()
{
    cerr << "make_earsdb input_dblistfile outdb" <<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    if(argc!=3) usage();
    ios::sync_with_stdio();
    FILE *fp;
    fp=fopen(argv[1],"r");
    if(fp==NULL)
    {
        cerr << "Open failed on file "<<argv[1]<<endl;
        usage();
    }
    string outdb(argv[2]);
    DatascopeHandle dbo(outdb,false);
    SEEDSiteTable stbl;
    MetadataList mdl;
    AttributeMap am("css3.0");
    EventCatalog origins(dbo,mdl,am);
    int count(0);
    char line[512];
    double t0min,t0max;  // valid range for origin times 
    t0min=epoch(1990001);
    t0max=epoch(2020001);
    /* Enclosing this entire loop in try block.  Not the best error
       handling scheme, but simplifies this. */
    try{
    while(fscanf(fp,"%s",line)!=EOF)
    {
        string dbname(line);
        Metadata originmd;  // empty metadata object
        cout << "Attempting to process leaf database with name "<<dbname<<endl;
        DatascopeHandle dbh(line,true);
        dbh.lookup("origin");
        dbh.rewind(); 
        /* Allow processing data organized by station or event.  For
           station version may have multiple origins per db while
           for events there is only one.  Will work either way because
           origins object will automatically handle duplicates. */
        long number_origins=dbh.number_tuples();
        for(long irec=0;irec<number_origins;++irec,++dbh)
        {
            double lat,lon,depth,t0;
            try {
                lat=dbh.get_double("lat");
                lon=dbh.get_double("lon");
                depth=dbh.get_double("depth");
                t0=dbh.get_double("time");
            }catch (SeisppError& serr)
            {
                cerr << "Error reading origin table for db "
                    << dbname<<endl<<"Error throw:"<<endl;
                serr.log_error();
                cerr << "Fix this db and try again"<<endl;
                exit(-1);
            }
            if( (lat>=-90.0)&&(lat<=90.0)
                    &&(lon>=-180.0)&&(lon<=360.0)
                    &&(t0>t0min) && (t0<t0max))
            {
                string method("tttaup"),model("iasp91");
                Hypocenter h(rad(lat),rad(lon),depth,t0,method,model);
                origins.add(h,originmd);
            }
            else
                cerr << "Bad hypocenter read: lat="<<lat
                    <<" lon="<<lon<<endl
                    <<"This row of db will be ignored"<<endl;
        }
            /* this was the older code before June 26 rewrite 
            dbh.lookup("site");
            This is the revision */
            dbh.lookup("wfdisc");
            dbh.natural_join("site");
            dbh.subset("chan=~/.*R/");
            double elev;
            long nrows=dbh.number_tuples();
            dbh.rewind();
            for(long i=0;i<nrows;++i,++dbh)
            {
                string net;
                /* Older code extracted net name from the dbname like this
                net=extract_net_string(dbname);
                Now we extract dir from the wfdisc row and crack that string */
                string wfdir=dbh.get_string("dir");
                net=extract_net_string(wfdir);
                //cout << "dir file="<<wfdir <<" extracted net="<<net<<endl;

                try {
                    SEEDStaGeom teststa(dbh.db,net);
                    if(coords_bad(teststa)) 
                    {
                        cout << "Warning bad coordinates found for net "
                            << net << " sta=" << teststa.sta<<endl
                            << "Found lat="<<teststa.lat<<" lon="<<teststa.lon<<endl
                            << "Dropped from site geometry"
                            <<endl;
                        continue;
                    }
                    //cout << teststa<<endl;
                    int retcode=stbl.add(teststa);
                    if(retcode<0)
                    {
                        SEEDStaGeom dupsta=stbl.last_duplicate();
                        if(teststa==dupsta)
                        {
                        cout << "Found duplicate net:sta with very different locations"<<endl;
                        cout << "Current:  "<<endl<<teststa<<endl;
                        cout << "Duplicate stored previously:  "
                            <<endl<<dupsta<<endl;
                        cout << "Replace or not?"<<endl
                            <<" Enter y to use the new (current). "
                            <<"Any other key will retain previous table values:";
                        string ques;
                        cin>>ques;
                        if(ques=="y")
                            stbl.replace(teststa);
                        }
                    }
                    else if(retcode>0 && SEISPP_verbose)
                    {
                        cout << "Found duplicate net:sta with a minor location difference"<<endl;
                        cout << "Current:  "<<teststa<<endl;
                        SEEDStaGeom dupsta=stbl.last_duplicate();
                        cout << "Stored:  "<<dupsta<<endl;
                    }
                }catch(SeisppError& serr)
                {
                    serr.log_error();
                    cerr << "Fix problem and try again"<<endl;
                    exit(-1);
                }
                dbh.close();
            }
        }
        /* Now we save origin, site, and snetsta tables in output db */
        int norigins=save_origin(dbo,origins);
        cout << "Saved "<<norigins<<" records in event->origin table"<<endl;
        /* save_site returns this container that defines renaming of stations
           required for duplicate names*/
        multimap<string,string> dups;  
        dups=save_site(dbo,stbl);
        cout << "Saved site, sitechan, and snetsta tables"<<endl
            << "Number of stations requiring name changes="<<dups.size()<<endl;
        /* Now build the wfdisc from all the pieces.  Change dir fields and
           (more importantly) the sta name when necessary */
        rewind(fp);
        string basedirname,dirname;
        dbo.lookup("wfdisc");
        while(fscanf(fp,"%s",line)!=EOF)
        {
            DatascopeHandle dbh(line,true);
            dbh.lookup("wfdisc");
            long nrows=dbh.number_tuples();
            /* This is is one level up from the final dir field.
            Peculiarity of the way these are assembled based on the
            ears dir structure. */
            basedirname=extract_directory_name(line);
            Dbptr db;
            db=dbo.db;
            char wfdrow[512];
            dbh.rewind();
            //DEBUG
            cout << "Attempting to write "<<nrows<<" rows of wfdisc from db "
                <<line<<endl;
            for(long i=0;i<nrows;++i,++dbh)
            {
                char chsta[10];
                char chdir[100];
                dbget(dbh.db,wfdrow);
                dbgetv(dbh.db,0,"sta",chsta,"dir",chdir,NULL);
                string sta(chsta);
                // previous file organization used this
                //dirname=basedirname+"/"+sta;
                /* New adds a net.sta naming in dir */
                string wfdir(chdir);
                dirname=basedirname+"/"+wfdir;
                const int MAXDIRSIZE(64);  // from css3.0 wfdisc 
                if(dirname.size()>MAXDIRSIZE) 
                {
                    cerr << "Error in writing wfdisc record "<<i<<endl
                        << "Directory name="<<dirname<<endl
                        << "Exceeds maximum length allowed in wfdisc of "<<MAXDIRSIZE<<endl
                        << "Edit directory name list to a relative path"<<endl
                        << "Cannot continue"<<endl;
                    exit(-4);
                }
                long rec=dbadd(dbo.db,wfdrow);
                if(rec==dbINVALID)
                {
                    cerr << "dbadd failed for row "<<i<<" of db "
                        <<line<<endl;
                }
                else
                {
                    db.record=rec;
                    dbputv(db,0,"dir",dirname.c_str(),NULL);
                    if(dups.count(sta)>1)
                    {
                        /* as above we need to change the 
                           way this procedure is called besause
                           naming changed 
                        string net=extract_net_string(line);
                        */
                        string net=extract_net_string(wfdir);
                        string seedsta=seedname(net,sta);
                        dbputv(db,0,"sta",seedsta.c_str(),NULL);
                    }
                }
            }
            dbh.close();
        }
    }catch(SeisppError& serr)
    {
        serr.log_error();
        cerr << "Fix problem and try again"<<endl;
        exit(-1);
    }
}
