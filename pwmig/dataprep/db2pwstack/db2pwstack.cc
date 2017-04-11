#include <stdio.h>
#include <string>
#include <list>
#include "seispp.h"
#include "TimeWindow.h"
#include "Metadata.h"
#include "ensemble.h"
#include "dbpp.h"
#include "pwstack_reader.h"
using namespace std;
using namespace SEISPP;
const string svmrowkey("SlownessVectorMatrixRows");
const string svmcolkey("SlownessVectorMatrixColumns");
int number_live_members(ThreeComponentEnsemble& gather)
{
    int live_count;
    vector<ThreeComponentSeismogram>::iterator gptr;
    for(gptr=gather.member.begin(),live_count=0;
		gptr!=gather.member.end();++gptr)
    {
	if(gptr->live)++live_count;
    }
    return live_count;
}
void LoadGatherHeader(ThreeComponentEnsemble& gather,PwstackGatherHeader& gh,
		int sequence_number)
{
    try {
	gh.sequence_number=sequence_number;
	int gather_size=gather.member.size();
        gh.number_members=number_live_members(gather);
        long evid=gather.get_long("evid");
        if(SEISPP_verbose) 
	{
            cout << "pwstack file writer: evid="<<evid<<" has "
                <<gh.number_members << " live three-component seismograms"<<endl;
	    if(gather_size!=gh.number_members)
		cout << (gather_size - gh.number_members) << " gather members "
			<< "are marked dead and will be deleted" <<endl;
	}
        gh.evid=evid;
        /* Warning - careful units are consistent with expectation of pwstack*/
        gh.lat=gather.get_double("origin.lat");
        gh.lon=gather.get_double("origin.lon");
        gh.depth=gather.get_double("origin.depth");
        gh.origin_time=gather.get_double("origin.time");
        /* These are slowness grid dimension written to a separate file */
        gh.svmrows=gather.get_int(svmrowkey);
        gh.svmcolumns=gather.get_int(svmcolkey);
        
    }catch(SeisppError& serr)
    {
        cerr << "LoadGatherHeader procedure:  Probably a problem with ensemble Metadata"
            << endl<<"Error caught:"<<endl;
        serr.log_error();
        cerr << "Contents of ensemble metadata"<<endl
            << dynamic_cast<Metadata&>(gather)<<endl;
        cerr << "Fatal Error - cannot recover"<<endl;
        exit(-1);
    }
}
void LoadTraceHeader(ThreeComponentSeismogram& d,PwstackTraceHeader& th)
{
    try {
        string sta=d.get_string("sta");
        strncpy(th.sta,sta.c_str(),8);
        th.time=d.t0;
        th.endtime=d.endtime();
        th.nsamp=d.ns;
        th.samprate=1.0/d.dt;
        th.atime=d.get_double("arrival.time");
        th.lat=d.get_double("site.lat");
        th.lon=d.get_double("site.lon");
        th.elev=d.get_double("site.elev");
    }catch(...){throw;};
}
/* This procedure creates a Hypocenter object from the
   ensemble Metadata using frozen keys */
Hypocenter LoadHypocenter(ThreeComponentEnsemble& d)
{
    try{
        double lon,lat,depth,time;
        lon=d.get_double("origin.lon");
        lat=d.get_double("origin.lat");
        depth=d.get_double("origin.depth");
        time=d.get_double("origin.time");
        /* here we have to convert these to radians */
        lon=rad(lon);
        lat=rad(lat);
        /* for now the method/model is frozen */
        Hypocenter h(lat,lon,depth,time,string("tttaup"),string("iasp91"));
        return(h);
    }catch(...){throw;};
}

/* This procedure loads slowness vectors ux,uy in the field
   variable of GCLvectorfield ug.  Travel times are computed
   with the hypocenter object h.  
*/
void BuildSlownessGrid(Hypocenter h, GCLvectorfield& ug)
{
    try {
        int i,j;
        double lat,lon;
        double elev(0.0);   // forced to sea level datum  - minor detail`
        for(i=0;i<ug.n1;++i)
            for(j=0;j<ug.n2;++j)
            {
                double lat,lon;
                lat=ug.lat(i,j);
                lon=ug.lon(i,j);
                /* Warning - for now this is frozen as P phase.  If this
                   is ever adapted to S to P data this needs to change*/
                SlownessVector u=h.pslow(lat,lon,elev);
                ug.val[i][j][0]=u.ux;
                ug.val[i][j][1]=u.uy;
//DEBUG
/*
cout << "source lat,lon="
<<deg(h.lat)<<", "<<deg(h.lon)
<< " psta lat,lon="
<< deg(lat)<<", "<<deg(lon)
<< " Back azimuth="<< deg(u.baz())<<endl;
*/
            }
    }catch(...){throw;};
}
/* This procedure takes slowness vector data stored in the 
   GCLvectorfield object u and puts the data in a linear buffer
   ubuf for writing.  It is blindly assumed ubuf has a size
   consistent with the dimensions of u.

Intentionally has no error handling. */
void LoadSlownessBuffer(GCLvectorfield& ug, double *ubuf)
{
        int i,j,k;
        k=0;
        for(i=0;i<ug.n1;++i)
            for(j=0;j<ug.n2;++j)
            {
                ubuf[k]=ug.val[i][j][0];
                ++k;
                ubuf[k]=ug.val[i][j][1];
                ++k;
            }
}
int write_ensemble(ThreeComponentEnsemble& g,FILE *fp)
{
    const string base_error("Error in write_ensemble procedure: ");
    try {
	if(SEISPP_verbose)
		cout << "Processing ensemble with the following attributes:"<<endl
			<< dynamic_cast<Metadata&>(g)<<endl;
        PwstackTraceHeader th;
        int n=g.member.size();
        vector<ThreeComponentSeismogram>::iterator gptr;
        int count;
        for(gptr=g.member.begin(),count=0;gptr!=g.member.end();++gptr,++count)
        {
	    if(gptr->live)
	    {
                LoadTraceHeader(*gptr,th);
                if(fwrite(&th,sizeof(PwstackTraceHeader),1,fp)!=1)
                    throw SeisppError(base_error
                            + "fwrite error writing trace header for station "
                            + th.sta);
                double *dptr=gptr->u.get_address(0,0);
                int nstotal=3*(gptr->ns);
                if(fwrite(dptr,sizeof(double),nstotal,fp)!=nstotal)
                    throw SeisppError(base_error
                            + "fwrite error writing sample data for station "
                            + th.sta);
	   }
           else if(SEISPP_verbose)
	   {
		/* Note this logic always skips data marked dead but only writes a message
		if verbose is set */
		cout << "Warning:   member number "<<count<<" of ensemble was marked dead"<<endl
		  << "Probably database problem.   Metadata for this seismogram follow:"<<endl
		  << dynamic_cast<Metadata&>(*gptr)<<endl;
           }
        }
        return(count);
    } catch(...){throw;};
}

            
void usage()
{
    cerr << "db2pwstack db outfile [-v -pf pffile]"<<endl;
    exit(-1);
}
void write_directory(long *ids,long *foffs, int nevents,FILE *fp)
{
    const string base_error("write_error procedure:  ");
    /* Probably should do a seek to end of file, but not necessary 
       with current logic */
    long diroffset=ftell(fp);
    if(fwrite(ids,sizeof(long),nevents,fp)!=nevents)
        throw(base_error + "fwrite failed writing id vector");
    if(fwrite(foffs,sizeof(long),nevents,fp)!=nevents)
        throw(base_error + "fwrite failed writing foff vector");
    rewind(fp);
    if(fwrite(&diroffset,sizeof(long),1,fp)!=1)
        throw(base_error + "fwrite failed writing director offset address");
    // Converted to long to avoid mixed types for this pair of values
    // written at top of the file
    long nevout=(long)nevents;
    if(fwrite(&nevout,sizeof(long),1,fp)!=1) 
        throw(base_error + "fwrite failed writing nevents value");
}

bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    ios::sync_with_stdio();
    if(argc<3) usage();
    string dbname(argv[1]);
    string outfile(argv[2]);
    string pffile("db2pwstack");
    int i;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pffile=string(argv[i]);
        }
        else if(sarg=="-v")
            SEISPP_verbose=true;
        else
            usage();
    }
    Pf *pf;
    if(pfread(const_cast<char *>(pffile.c_str()),&pf))
    {
        cerr << "pfread failed for pffile="<<pffile<<endl;
        usage();
    }
    FILE *fp;
    fp=fopen(outfile.c_str(),"w");
    if(fp==NULL) 
    {
        cerr << "fopen failed trying to open output file="<<outfile<<endl;
        exit(-1);
    }
    try{
        DatascopeHandle dbh(dbname, true);
        Metadata control(pf);
        string arrival_phase=control.get_string("phase_key");
	int minimum_gather_size=control.get_int("minimum_gather_size");
	string dbviewmode(control.get_string("database_view_mode"));
	if(dbviewmode=="dbprocess")
        	dbh=DatascopeHandle(dbh,pf,string("dbprocess_commands"));
	else if(dbviewmode=="use_wfdisc")
	{
		dbh.lookup("arrival");
		list<string> j1,j2;
		j1.push_back("sta");
		j1.push_back("wfdisc.time::wfdisc.endtime");
		j2.push_back("sta");
		j2.push_back("arrival.time");
		dbh.leftjoin("wfdisc",j1,j2);
		dbh.natural_join("assoc");
                dbh.subset("phase=~/"+arrival_phase+"/");
		dbh.natural_join("origin");
		dbh.natural_join("event");
		dbh.subset("orid==prefor");
		j1.clear();
		j2.clear();
		j1.push_back("sta");
		j1.push_back("chan");
		j1.push_back("arrival.time");
		j2.push_back("sta");
		j2.push_back("chan");
		j2.push_back("ondate::offdate");
		dbh.join("sitechan",j1,j2);
		dbh.natural_join("site");
		if(SEISPP_verbose) cout << "working view size="
			<<dbh.number_tuples()<<endl;
		list<string> sortkeys;
		sortkeys.push_back("evid");
		sortkeys.push_back("sta");
		sortkeys.push_back("chan");
		dbh.sort(sortkeys);
		list<string> gkey;
		gkey.push_back("evid");
		gkey.push_back("sta");
		dbh.group(gkey);
	}
	else if(dbviewmode=="use_wfprocess")
	{
		dbh.lookup("event");
		dbh.natural_join("origin");
		dbh.subset("orid==prefor");
		dbh.natural_join("assoc");
		dbh.natural_join("arrival");
                dbh.subset("phase=~/"+arrival_phase+"/");
                if(SEISPP_verbose) cout << "Catalog view size="<<dbh.number_tuples()<<endl;
		DatascopeHandle ljhandle(dbh);
		ljhandle.lookup("wfprocess");
		ljhandle.natural_join("sclink");
		ljhandle.natural_join("evlink");
                if(SEISPP_verbose)cout << "Left join table size="<<ljhandle.number_tuples()<<endl;
		list<string> jk,sjk;
		jk.push_back("evid");
		jk.push_back("sta");
		dbh.join(ljhandle,jk,jk);
                if(SEISPP_verbose)cout << "Working table size="<<dbh.number_tuples()<<endl;
                sjk.push_back("sta");
		dbh.join(string("site"),sjk,sjk);
                if(SEISPP_verbose)cout << "After site join size="<<dbh.number_tuples()<<endl;
		list<string> sortkeys;
		sortkeys.push_back("evid");
		sortkeys.push_back("sta");
		dbh.sort(sortkeys);
	}
	else
	{
		cerr << "Illegal option for parameter database_view_mode="<<dbviewmode;
		exit(-1);
	}
        MetadataList station_mdl=pfget_mdlist(pf,"station_metadata");
        MetadataList ensemble_mdl=pfget_mdlist(pf,"ensemble_metadata");
        AttributeMap InputAM("css3.0");
        double tsfull, tefull;
        tsfull = control.get_double("data_time_window_start");
        tefull = control.get_double("data_time_window_end");
        TimeWindow data_window(tsfull,tefull);
        string gridname=control.get_string("PseudostationGridName");
        bool psg_use_file=control.get_bool("use_file_for_pseudostationgrid");
        GCLgrid *g;
        if(psg_use_file)
            g=new GCLgrid(gridname);
        else
            g=new GCLgrid(dbh,gridname);
        /* Slowness vector data stored here */
        GCLvectorfield u(*g,2);
        /* this is used as a write buffer to simplify writing of svm data */
        int nu=2*u.n1*u.n2;  // 2 because of ux,uy pairs for each grid point
        double *ubuf=new double[nu];  // not freed in this code - warning
        delete g;
	cout << "Processing begins on database " 
		<<  dbname << endl
		<<"Number of rows in working database view== "<<dbh.number_tuples() <<endl;
        list<string> group_keys;
        group_keys.push_back("evid");
        dbh.group(group_keys);
        dbh.rewind();
        int nevents=dbh.number_tuples();
	cout << "This run will process data from "<<nevents<<" events."<<endl;
        PwstackGatherHeader gh;
        long *ids=new long[nevents];
        long *foffs=new long[nevents];
        long diroffset(0);
        /* This is a dummy write for now - filled in as last step*/
        fwrite(&diroffset,sizeof(long),1,fp);
        fwrite(&diroffset,sizeof(long),1,fp);
        int rec;
	int number_events_processed=dbh.number_tuples();
	int number_events_saved(0);
        for(rec=0,dbh.rewind();rec<dbh.number_tuples();++rec,++dbh)
        {
            /* This requires segmented data */ 
        	ThreeComponentEnsemble *din = new ThreeComponentEnsemble(dynamic_cast<DatabaseHandle&>(dbh),
	                                         station_mdl, ensemble_mdl,InputAM);
	        auto_ptr<ThreeComponentEnsemble> ensemble = ArrivalTimeReference(*din,
	                "arrival.time",data_window);
	        delete din;
            long evid=ensemble->get_long("evid");
	    int gather_size=ensemble->member.size();
	    cout << "Event id="<<evid<<endl
		<< "source latitude="<<ensemble->get_double("origin.lat")
		<< " longitude="<<ensemble->get_double("origin.lon")
		<< " depth="<< ensemble->get_double("origin.depth")
		<< " origin time="
			<< strtime(ensemble->get_double("origin.time"))
		<<endl;
	    if(gather_size<minimum_gather_size)
	    {
		cout << "Warning:  number of seismograms="
			<< gather_size<<" is below threshold of "
			<< minimum_gather_size<<endl
			<< "Data for this event not written to output file"
			<<endl;
	    }
	    else
	    {
                Hypocenter h(LoadHypocenter(*ensemble));
                BuildSlownessGrid(h,u);
                LoadSlownessBuffer(u,ubuf);
                /* Note row and column are confusing concepts here because
                   of C and fortran difference.  row is index 1 */
                ensemble->put(svmrowkey,u.n1);
                ensemble->put(svmcolkey,u.n2);
                LoadGatherHeader(*ensemble,gh,number_events_saved);
                long pos=ftell(fp);
                ids[number_events_saved]=evid;
                foffs[number_events_saved]=pos;
                if(fwrite(&gh,sizeof(PwstackGatherHeader),1,fp)!=1)
                {
                    cerr << "fwrite failed writing ensemble header event for evid="
                        <<evid<<endl;
                    exit(-1);
                }
                if(fwrite(ubuf,sizeof(double),nu,fp)!=nu)
                {
                    cerr << "fwrite failed writing slowness vector matrix section "
                        << "for evid="<<evid<<endl;
                    exit(-1);
                }
                int nseis=write_ensemble(*ensemble,fp);
		cout << "Wrote "<<nseis<<" seismograms of expected "
                    << gather_size<<endl;
		++number_events_saved;
	    }
        }
	/* Note we use number_events_saved because above loop
	can and usually will drop some gathers. write_directory
	and the file structure handles this correctly */
        write_directory(ids,foffs,number_events_saved,fp);
        delete [] ids;
        delete [] foffs;
    }catch(std::exception& excp)
    {
        cerr << excp.what()<<endl;
    }
    catch(SeisppError& serr)
    {
        serr.log_error();
    }
}
