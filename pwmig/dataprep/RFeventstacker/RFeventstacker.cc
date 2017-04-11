#include <exception>
#include <algorithm>
#include "seispp.h"
#include "stack.h"
#include "Hypocenter.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
#include "Metadata.h"
#include "dbpp.h"
#include "VectorStatistics.h"
#include "resample.h"
#include "SimpleWavelets.h"
using namespace std;
using namespace SEISPP;
/* This procedure will eventually be in libseispp under resample.h but for 
   testing need the prototype here.*/
auto_ptr<ThreeComponentEnsemble> Resample(ThreeComponentEnsemble& d, 
        ResamplingDefinitions& rd, double dtout, bool trim);
/* This is a simplified version of a template in XcorProcessingEngine.  Simpler
is cleaner in this case*/
template<class T> struct less_stackweight
                : public binary_function<T,T,bool> {
	bool operator()(T x, T y)
        {
		string keyword=SEISPP::stack_weight_keyword;
               double valx,valy;
               try {
			valx=x.get_double(keyword);
		}catch(...){return true;};
		try{
                        valy=y.get_double(keyword);
                }catch(...) {return false;};

		if(valx<valy)
	                return true;
	        else
	                return false;
	 }
};



/* This procedure does not actually seem to be necessary.  Think we
can just build a view containing the hypocentroid and we should be fine.  
Will retain this, however, until I'm sure that is so */

/***************************
EventCatalog LoadHypocentroids(DatascopeHandle dbh,
    string gridname,
	string firstotime, double otoffset,ttmethod,ttmodel)
{
    const string base_error("LoadHypocentroid:  ");
    try {
	double ot0=str2epoch(firstotime.c_str*());
	double hclat,hclon,hcdepth,otimei(ot0);
	MetadataList mdl;
	const Metadata_typedef mdaux={string("gridid"),MDint};
	mdl.push_back(mdaux);
	EventCatalog result(mdl);
	dbh.lookup("hypocentroid");
	dbh.subset(string("gridname=~/")+gridname+"/";
	dbh.rewind();
	int nrec=dbh.number_tuples();
	if(nrec<=0) throw SeisppError(base_error
			+ "hypocentroid table has no data for table"
			+gridname);
	int i,gridid,nass;
	Metadata aux;
	for(i=0;i<nrec;++i,++dbh,otime+=otoffset)
	{
		lat=dbh.get_double("hclat");
		lon=dbh.get_double("hclon");
		depth=dbh.get_double("hcdepth");
		gridid=dbh.get_int("gridid");
		nass=dbh.get_int("nass");
		aux.put("gridid",gridid);
		aux.put("nass",nass);
		Hypocenter h(lat,lon,depth,otime,ttmethod,ttmodel);	
		if(!result.add(h,aux))
		{
			throw SeisppError(base_error
				+ "hypocentroid has duplicates.  fix db");
		}
	}
	// loop over db, load hypo data, set otime, then call "add" method
	// of EventCatalog -- one for each hypo
    } catch (...) {throw;}
}
**********************/

Hypocenter MakeHypocentroid(ThreeComponentEnsemble& d,double otime,
	string ttmethod, string ttmodel)
{
    try {
    	double lat,lon,depth;
	lat=d.get_double("hclat");
	lon=d.get_double("hclon");
	depth=d.get_double("hcdepth");
	Hypocenter h(rad(lat),rad(lon),depth,otime,ttmethod,ttmodel);
	return(h);
    } catch (...){throw;};
}

DatascopeHandle BuildWaveformView_wfprocess(DatascopeHandle& dbhin,string phase)
{
	const string base_error("RFEventStacker::BuildWaveformView_wfprocess(Fatal Error):  ");
	try {
		DatascopeHandle dbh(dbhin);
		dbh.lookup("hypocentroid");
		dbh.natural_join("cluster");
		dbh.natural_join("event");
		dbh.natural_join("origin");
		dbh.subset("orid==prefor");
		dbh.natural_join("assoc");
		dbh.natural_join("arrival");
		string phase_subset;
		phase_subset="phase=~/"+phase+"/";
		dbh.subset(phase_subset);
		if(dbh.number_tuples()<=0) throw SeisppError(base_error
			+ "hypocentroid->cluster->catalogdata"
			+ "view has no data");
		if(SEISPP_verbose) cout << "hypocentroid->cluster->catalogdata "
			<< "view size="<<dbh.number_tuples()<<endl;
		DatascopeHandle ljhandle(dbh);
		ljhandle.lookup("wfprocess");
		ljhandle.natural_join("sclink");
		ljhandle.natural_join("evlink");
		if(ljhandle.number_tuples()<=0) throw SeisppError(base_error
			+ "wfprocess->sclink->evlink join has no data");
		if(SEISPP_verbose) cout << "waveform view size="
					<< ljhandle.number_tuples()<<endl;
		list<string> jk;
		jk.push_back("evid");
		jk.push_back("sta");
		dbh.join(ljhandle,jk,jk);
		list<string> sitejoinkeys;
		sitejoinkeys.push_back("sta");
		sitejoinkeys.push_back("ondate::offdate");
		list<string> viewjkeys;
		viewjkeys.push_back("sta");
		viewjkeys.push_back("wfprocess.time");
		dbh.join("site",viewjkeys,sitejoinkeys);
		if(dbh.number_tuples()<=0) throw SeisppError(base_error
			+ "join of left and right db tables has no data");
		if(SEISPP_verbose) cout << "Full working waveform view size="
					<< dbh.number_tuples()<<endl;
		dbh.subset("datatype=~/3c/ || datatype=~/c3/");
		if(dbh.number_tuples()<=0) throw SeisppError(base_error
			+ "subset for 3c bundle data failed.  Check wfprocess table");
		if(SEISPP_verbose) cout << "View size of 3c data="
					<< dbh.number_tuples()<<endl;
		list<string> sortkeys;
		sortkeys.push_back("gridid");
		sortkeys.push_back("sta");
		sortkeys.push_back("evid");
		dbh.sort(sortkeys);
		list<string> group_keys;
		group_keys.push_back("gridid");
		group_keys.push_back("sta");
		dbh.group(group_keys);
		if(SEISPP_verbose) cout << "Number of waveform ensembles "
					<< "(Number of stacks to compute) = "
					<< dbh.number_tuples()<<endl;
		dbh.rewind();
		return(dbh);
	} catch (...) {throw;};
}
DatascopeHandle BuildWaveformView_wfdisc(DatascopeHandle& dbhin,string phase)
{
    try{
        const string base_error("RFEventStacker::BuildWaveformView_wfdisc(Fatal Error):  ");
	DatascopeHandle dbh(dbhin);
	dbh.lookup("hypocentroid");
	dbh.natural_join("cluster");
	dbh.natural_join("event");
	dbh.natural_join("origin");
	dbh.subset("orid==prefor");
	dbh.natural_join("assoc");
	dbh.natural_join("arrival");
	string phase_subset;
	phase_subset="phase=~/"+phase+"/";
	dbh.subset(phase_subset);
	if(dbh.number_tuples()<=0) throw SeisppError(base_error
		+ "hypocentroid->cluster->catalogdata"
		+ "view has no data");
	if(SEISPP_verbose) cout << "hypocentroid->cluster->catalogdata "
		<< "view size="<<dbh.number_tuples()<<endl;
        list<string> j1,j2;
        j1.push_back("sta");
        j1.push_back("wfdisc.time::wfdisc.endtime");
        j2.push_back("sta");
        j2.push_back("arrival.time");
        dbh.leftjoin(string("wfdisc"),j1,j2);
        if(SEISPP_verbose) cout << "size after wfdisc join="
            << dbh.number_tuples()<<endl;
        dbh.natural_join(string("site"));
        /* Natural join of sitechan is problematic because it uses
           chanid.   This has been problematic since dogmatic 
           requirement of sensor table to run dbfixchanids.   
           Consequently have to do this with explicit keys */
        j1.clear();
        j2.clear();
        j1.push_back("sta");
        j1.push_back("chan");
        j2=j1;
        j1.push_back("arrival.time");
        j2.push_back("ondate::offdate");
        dbh.join(string("sitechan"),j1,j2);
        if(SEISPP_verbose) cout << "working view size="
            << dbh.number_tuples()<<endl;
        list<string> sortkeys;
        sortkeys.push_back("gridid");
        sortkeys.push_back("sta");
        sortkeys.push_back("evid");
        sortkeys.push_back("chan");
        dbh.sort(sortkeys);
        list<string> gkey;
        gkey.push_back("gridid");
        gkey.push_back("sta");
        gkey.push_back("evid");
        dbh.group(gkey);
        if(SEISPP_verbose) cout << "Number of three-component"
            <<" bundles in working view="<<dbh.number_tuples()<<endl;
        gkey.clear();
        gkey.push_back("gridid");
        gkey.push_back("sta");
        dbh.group(gkey);
        dbh.rewind();
        if(SEISPP_verbose) cout << "Number of ensembles in working view="
            <<dbh.number_tuples()<<endl;
        return(dbh);
    }catch(...){throw;};
}
/* exit with a message if a table has any existing content.
Required for this program.  Bad bad form, but moves us forward.
Ultimately this will be replace by output to a file or an HPC db 
that is not relational */
void exit_if_exists(DatascopeHandle& dbh,char *table)
{
	int nrec=dbh.number_tuples();
	if(nrec>0)
	{
	    cerr << "RFeventstacker:  output database table "
		<< table << " is not empty"<<endl
		<< "This program requires an empty output db"<<endl
		<< "Change output name or remove old tables"<<endl;
	    exit(-1);
	}
}
/* necessary for current 3c seismogram constructor.  This is very 
crude and assuems components passed in d are cardinal components 
in e,n,z order */
void load_hang_vang(vector<TimeSeries>& d)
{
	d[0].put("hang",90.0);
	d[0].put("vang",90.0);
	d[1].put("hang",0.0);
	d[1].put("vang",90.0);
	d[2].put("hang",0.0);
	d[2].put("vang",0.0);
}
/* These two procedures are special for this program.  They write a minimal
amount of information ot arrival and assoc. */
long int save_arrival(DatascopeHandle& dbh, ThreeComponentSeismogram& d,
		string arrivalchan, string phase)
{
    try{
	double atime=d.get_double("arrival.time");
	string sta=d.get_string("sta");
	long jday=yearday(atime);
	dbh.append();
	dbh.put("sta",sta);
	dbh.put("time",atime);
	long arid;
	arid = dbnextid(dbh.db,"arid");
	dbh.put("arid",arid);
	dbh.put("jdate",jday);
	dbh.put("chan",arrivalchan);
	dbh.put("iphase",phase);
	return(arid);
    }catch (...) {throw;};
}
void save_assoc(DatascopeHandle& dbh,ThreeComponentSeismogram& d, 
		int orid, int arid, string phase, Hypocenter& h)
{
    try{
	string sta=d.get_string("sta");
	double stalat,stalon,staelev;
	stalat=d.get_double("site.lat");
	stalon=d.get_double("site.lon");
	staelev=d.get_double("site.elev");
	double delta=h.distance(stalat,stalon);
	double seaz=h.seaz(stalat,stalon);
	double esaz=h.esaz(stalat,stalon);
	string timedef("d");
	dbh.append();
	dbh.put("arid",arid);
	dbh.put("orid",orid);
	dbh.put("sta",sta);
	dbh.put("phase",phase);
	dbh.put("delta",delta);
	dbh.put("seaz",seaz);
	dbh.put("esaz",esaz);
	dbh.put("timedef",timedef);
	dbh.put("timeres",0.0);
    }catch (...) {throw;};
}
double coherence(double ssqr,double ssqd)
{
    // This is an ambiguous choice to avoid a NaN but
    // I make it 0 as a safer choice.  Avoids 0/0
    if(ssqr<=0.0 && ssqd <= 0.0)
        return(0.0);
    else if(ssqr>ssqd)
        return(0.0);
    else
        return(1.0-sqrt(ssqr)/sqrt(ssqd));
}

typedef struct {
	double coherence;
	double semblance;
} StackStatistics;
StackStatistics ComputeCoherence(TimeSeriesEnsemble& d1in, TimeSeriesEnsemble& d2in, TimeSeriesEnsemble& d3in, 
	Stack& s1, Stack& s2, Stack& s3, TimeWindow twin, bool iz)
{
	/* We can assume d1,d2, and d3 are all the same size here. Don't blindly use this elsewhere without noting this*/
	int n=d1in.member.size();
	double ssr,ssd,resid;
	double ssqstack;
	int i,j,ns;
	for(ssqstack=0.0,i=0;i<s1.stack.ns;++i)
	{
		ssqstack+=(s1.stack.s[i]*s1.stack.s[i]);
		ssqstack+=(s2.stack.s[i]*s2.stack.s[i]);
		if(!iz)ssqstack+=(s3.stack.s[i]*s3.stack.s[i]);
	}
	for(ssr=0.0,ssd=0.0,i=0;i<n;++i)
	{
		if(!d1in.member[i].live) continue;
		TimeSeries d1=WindowData(d1in.member[i],twin);
		/* this is mostly for debug.  As written here this should not happend, but best leave this
		in as a sanity check.  Also intentionally only test d1 since d2 and d3 are totally parallel
		in structure */
		if(abs(d1.ns-s1.stack.ns)>1)
		{
			cerr << "Size mismatch between data ensemble and stack.  "<<endl
				<< "Stack number samples="<<s1.stack.ns
				<< " Windowed data size="<<d1.ns<<endl;
			exit(-1);
		}
		ns=min(s1.stack.ns,d1.ns);
		for(j=0;j<ns;++j)
		{
			ssd+=pow(d1.s[j],2);
			resid=d1.s[j] - s1.stack.s[j];
			ssr+=(resid*resid);
		}
		TimeSeries d2=WindowData(d2in.member[i],twin);
		ns=min(s2.stack.ns,d2.ns);
		for(j=0;j<ns;++j)
		{
			ssd+=pow(d2.s[j],2);
			resid=d2.s[j] - s2.stack.s[j];
			ssr+=(resid*resid);
		}
		if(!iz)
		{
			TimeSeries d3=WindowData(d3in.member[i],twin);
			ns=min(s3.stack.ns,d3.ns);
			for(j=0;j<ns;++j)
			{
				ssd+=pow(d3.s[j],2);
				resid=d3.s[j] - s3.stack.s[j];
				ssr+=(resid*resid);
			}
		}
	}
        /* Throw an exception if ssq of the data vector is zero.  Will get nan
           if this happens so is an exception.  Such data need to be dropped as it 
           indicates an empty data ensemble*/
        if(ssd<=0.0) throw SeisppError (string("ComputeCoherence:  ")
                + "ensemble in data window is all zero so cannot computer statistics");
	StackStatistics result;
	result.coherence=coherence(ssr,ssd);
	result.semblance=ssqstack*static_cast<double>(n)/ssd;
	return(result);
}
StackType GetStackType(Metadata& control)
{
    StackType result;
    try {
	string name=control.get_string("stack_type");
	if((name=="simple") || (name=="average") || (name=="straight_stack"))
	{
		result=BasicStack;
	}
	else if(name=="median")
	{
		result=MedianStack;
	}
	else if((name=="robust") || (name=="RobustSNR"))
	{
		result=RobustSNR;
	}
	else
	{
		cerr<<"WARNING:  Unrecognized tag for parameter stack_type="
			<< name<<endl
			<<"Defaulting to standard stack (mean)"<<endl;
		result=BasicStack;
	}
    } catch (MetadataGetError mderr)
    {
    	cerr << "Parameter stack_type is missing from input parameter file"
		<<endl
		<<"Defaulting to standard (mean) stack"<<endl;
	result=BasicStack;
    }
    return(result);
}
bool need_to_resample(ThreeComponentEnsemble& d,double targetdt)
{
    const double tolerance(0.001);  //allows soft equality test
    vector<ThreeComponentSeismogram>::iterator dptr;
    for(dptr=d.member.begin();dptr!=d.member.end();++dptr)
        if(fabs( ((dptr->dt)-targetdt)/targetdt)>tolerance) return true;
    return false;
}
TimeSeries build_filter_wavelet(Metadata& control)
{
    try {
        string filter_type=control.get_string("filter_type");
        int wavelet_length=control.get_int("wavelet_length");
        double width_parameter=control.get_double("wavelet_width_parameter");
        double dt=control.get_double("target_sample_interval");
        if(filter_type=="ricker")
        {
            return(ricker_wavelet(wavelet_length,dt,width_parameter,PEAK));
        }
        else if (filter_type=="gaussian")
        {
            return(gaussian_wavelet(wavelet_length,dt,width_parameter,AREA));
        }
        else
            throw SeisppError("build_filter_wavelet:   filter_type="
                    + filter_type + " not allowed\n" 
                    + "Must be either gaussian or ricker");
    }catch(...){throw;};
}
void usage()
{
	cerr << "RFeventstacker dbin dbout [-noplots -v -pf pfname]" << endl
		<< "dbout must be empty"<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
/* Major modification history:   
   Nov 14, 2012:   Originally pseudo events that are generated by this code from
   common receiver gathers of sources from common regions (telecluster) got orid and 
   evid values by simply counting.   In order to merge data from different generations 
   with conflicting problems, this was problematic.  To solve this I now make evid and
   orid be the same as the parent gridid defined by telecluster.  */
int main(int argc, char **argv)
{
	ios::sync_with_stdio();
	if(argc<3) usage();
	string pfin(argv[0]);
	string dbinname(argv[1]);
	string dboutname(argv[2]);
	bool plotdata(true);
        for(int i=3;i<argc;++i)
        {
                if(!strcmp(argv[i],"-v"))
		{
                        SEISPP_verbose=true;
		}
                else if(!strcmp(argv[i],"-pf"))
                {
                        ++i;
                        if(i>=argc) usage();
                        pfin = string(argv[i]);
                }
                else if(!strcmp(argv[i],"-noplots"))
		{
			plotdata=false;
		}
                else
                        usage();
        }
	Pf *pf;
        if(pfread(const_cast<char *>(pfin.c_str()),&pf)) 
	{
		cerr << "pfread error for pf file="<<pfin<<".pf"<<endl;
		exit(-1);
	}
	MetadataList mdlin=pfget_mdlist(pf,"Station_mdlist");
	MetadataList mdens=pfget_mdlist(pf,"Ensemble_mdlist");
	MetadataList mdwfprocess=pfget_mdlist(pf,"wfprocess_mdlist");
	MetadataList mdwfdisc=pfget_mdlist(pf,"wfdisc_mdlist");
	try {

                ResamplingDefinitions rd(pf);
		Metadata control(pf);
                /* The data will all be resampled to this sample interval */
                double dt0=control.get_double("target_sample_interval");
		/* This perhaps should be a parameter, but frozen for now*/
		const string atkey("arrival.time");
		double ts,te;
		ts=control.get_double("stack_starttime");
		te=control.get_double("stack_endtime");
		TimeWindow twin(ts,te);
		StackType stacktype=GetStackType(control);
		double rts,rte;
		rts=control.get_double("robust_window_starttime");
		rte=control.get_double("robust_window_endtime");
		TimeWindow rtwin(rts,rte);
		string phase=control.get_string("phase_for_alignment");
		string arrivalchan=control.get_string("arrival_chan");
		bool ignore_vertical=control.get_bool("ignore_vertical");
                /* These control an optional post stack convolution 
                   with a Ricker wavelet*/
                bool filter_stack=control.get_bool("filter_stack");
                TimeSeries filter_wavelet;
                if(filter_stack) filter_wavelet=build_filter_wavelet(control);

		/* All data will be written to this directory and to 
		one file set by dfile.  Intentional to improve performance
		in hpc systems*/
		string dir=control.get_string("output_dir");
		string dfile_base=control.get_string("output_dfile_base");
                string dfile;
                bool use_wfdisc=control.get_bool("use_wfdisc");
		DatascopeHandle dbh(string(dbinname),true);
                DatascopeHandle dbhwf(dbh);
                if(use_wfdisc)
		    dbhwf=BuildWaveformView_wfdisc(dbh,phase);
                else
		    dbhwf=BuildWaveformView_wfprocess(dbh,phase);
		if(dbhwf.number_tuples()<1)
		{
			cerr << "Waveform view has no data"<<endl
				<< "This join must be defined:  "
				<< "(wfprocess->sclink)<-"
				<< "hypocentroid->cluster->event->origin->assoc->arrival"
				<<endl;
			usage();
		}
		dbhwf.rewind();
		DatascopeHandle dbho(string(dboutname),false);
		DatascopeHandle dborigin(dbho);
		dborigin.lookup("origin");
		/* This procedure will exit the program if any of these
		tables have existing content */
		exit_if_exists(dborigin,"origin");
		DatascopeHandle dbevent(dbho);
		dbevent.lookup("event");
		exit_if_exists(dbevent,"event");
		DatascopeHandle dbassoc(dbho);
		dbassoc.lookup("assoc");
		exit_if_exists(dbassoc,"assoc");
		DatascopeHandle dbarrival(dbho);
		dbarrival.lookup("arrival");
		exit_if_exists(dbarrival,"arrival");
		DatascopeHandle dbsclink(dbho);
		dbsclink.lookup("sclink");
		exit_if_exists(dbsclink,"sclink");
		DatascopeHandle dbevlink(dbho);
		dbevlink.lookup("evlink");
		exit_if_exists(dbevlink,"evlink");
		DatascopeHandle dbwfprocess(dbho);
		dbwfprocess.lookup("wfprocess");
		exit_if_exists(dbwfprocess,"wfprocess");
		DatascopeHandle dbwfdisc(dbho);
		dbwfdisc.lookup("wfdisc");
		exit_if_exists(dbwfdisc,"wfdisc");
		DatascopeHandle dbstackstats(dbho);
		dbstackstats.lookup("stackstats");
		exit_if_exists(dbstackstats,"stackstats");

		const string schema("css3.0");
		AttributeMap am(schema);
		/* This builds a fake catalog of hypocentroids used to 
		guarantee proper associations */
		string firstotime,ttmethod,ttmodel;
		firstotime=control.get_string("first_fake_origin_time");
		double otoffset=control.get_double("fake_origin_time_offsets");
		ttmethod=control.get_string("ttmethod");
		ttmodel=control.get_string("ttmodel");
		double lat,lon,depth,otime;
		otime=str2epoch(const_cast<char *>(firstotime.c_str()));
		auto_ptr<TimeSeriesEnsemble> x1,x2,x3;  // hold components
		vector<string> chanmap;
		chanmap.push_back("E"); 
		chanmap.push_back("N"); 
		chanmap.push_back("Z"); 
		vector<TimeSeries> stack3c;
		stack3c.reserve(3);
		int nrec=dbhwf.number_tuples();
//DEBUG
cout << "Number of record to process="<<nrec<<endl;
		int evid,orid,arid;
		int record;
		int gridid,lastgridid;
		Hypocenter hcen;
		ThreeComponentEnsemble *pwdataraw;
		cout << "Station gridid fold coherence semblance"<<endl;
		for(record=0;record<nrec;++record,++dbhwf)
		{
			int fold;
                        stringstream ss;
                        /*
			pwdataraw=new ThreeComponentEnsemble(dynamic_cast<DatabaseHandle&>(dbhwf),
				mdlin,mdens,am);
                                */
                        auto_ptr<ThreeComponentEnsemble> pwdataraw =
                            auto_ptr<ThreeComponentEnsemble>(new
                                    ThreeComponentEnsemble(dbhwf,mdlin,mdens,am));
                        if(pwdataraw->member.size()<=0) 
                        {
                            cerr << "Ensemble for data for row "<<record<<" of database view is empty"<<endl
                                << "Likely missing waveform data or corrupted database"<<endl;
                            continue;
                        }
                        /*Scan for irregular sampling for efficiency.  Note always call
                          resampler with trim off assuming deconvolved data have no issue
                          with dc offset*/
                        if(need_to_resample(*pwdataraw,dt0))
                            pwdataraw=Resample(*pwdataraw,rd,dt0,false);
                        else
                        {
                            /* When we don't resample it is prudent to force dt of
                               all the data to the same exact value */
                            vector<ThreeComponentSeismogram>::iterator dptr;
                            for(dptr=pwdataraw->member.begin();
                                    dptr!=pwdataraw->member.end();++dptr) dptr->dt=dt0;
                        }
			gridid=pwdataraw->get_int("gridid");
			if(record==0) 
			{
				lastgridid=gridid;
				hcen=MakeHypocentroid(*pwdataraw,otime,
					ttmethod,ttmodel);
			}
                        /* The contents pushed to ss are used for evid */
                        ss<<lastgridid;
                        evid=lastgridid;
                        orid=lastgridid;
			if( (gridid!=lastgridid) || (record==(nrec-1)) )
			{
				dbevent.append();
				dbevent.put("evid",evid);
				dbevent.put("prefor",orid);
				dborigin.append();
				dborigin.put("lat",deg(hcen.lat));
				dborigin.put("lon",deg(hcen.lon));
				dborigin.put("depth",hcen.z);
				dborigin.put("time",hcen.time);
				dborigin.put("orid",orid);
				dborigin.put("evid",evid);
				/* Now we load the next hypo data */
				otime+=otoffset;
				hcen=MakeHypocentroid(*pwdataraw,otime,
					ttmethod,ttmodel);
				lastgridid=gridid;
			}
			string sta=pwdataraw->get_string("sta");
			auto_ptr<ThreeComponentEnsemble> 
			  pwdata = ArrivalTimeReference(*pwdataraw,atkey,twin);
			//delete pwdataraw;
			/* This appears to be necessary */
			for(int im=0;im<pwdata->member.size();++im)
				pwdata->member[im].put(moveout_keyword,0.0);
			x1=ExtractComponent(*pwdata,0);
			x2=ExtractComponent(*pwdata,1);
			x3=ExtractComponent(*pwdata,2);
			Stack *s1=NULL,*s2=NULL,*s3=NULL;
			try {
				/* This branching is not strictly necessary
				but is more efficient. */
				if(stacktype==BasicStack)
				{
					s1 = new Stack(*x1,twin);
					s2 = new Stack(*x2,twin);
					s3 = new Stack(*x3,twin);
				}
				else
				{
					s1=new Stack(*x1,twin,rtwin,stacktype);
					s2=new Stack(*x2,twin,rtwin,stacktype);
					s3=new Stack(*x3,twin,rtwin,stacktype);
				}
				fold=s1->fold;
                                StackStatistics stackstat;
				if(fold<=1)
				{
                                    stackstat.coherence=1.0;
                                    stackstat.semblance=1.0;
				}
				else
				{
                                    try {
				        stackstat=ComputeCoherence(*x1,*x2,*x3,
					        *s1,*s2,*s3,twin,ignore_vertical);
                                     } catch(SeisppError& serr)
                                     {
                                            cerr << "Error:  Station="<<sta<<" and gridid="<<gridid<<endl;
                                            serr.log_error();
                                            cerr << "Skipping these data."<<endl;
                                            continue;
                                     }
                                }
				cout << sta << " "
					<< gridid <<" "
					<< fold <<" "
					<< stackstat.coherence<<" "
					<< stackstat.semblance<<endl;
				dbstackstats.append();
				dbstackstats.put("gridid",gridid);
				dbstackstats.put("coherence",stackstat.coherence);
				dbstackstats.put("semblance",stackstat.semblance);
				dbstackstats.put("fold",fold);
				dbstackstats.put("sta",sta);
				dbstackstats.put("pchan","3c");
				dbstackstats.put("phase","P");
				stack3c.clear();
				stack3c.push_back(s1->stack);
				stack3c.push_back(s2->stack);
				stack3c.push_back(s3->stack);
				load_hang_vang(stack3c);
			}
			catch (SeisppError& serr)
			{
				cerr << "Stack constructor failed for sta="<<sta<<" and gridid="<<gridid<<endl
					<<"Stack returned this error message:"<<endl;
				serr.log_error();
				cerr << "Data for this station-gridid will be skipped"<<endl;
				continue;
			} 
			ThreeComponentSeismogram result(stack3c);
			if(s1!=NULL) delete s1;
			if(s2!=NULL) delete s2;
			if(s3!=NULL) delete s3;
                        /* This is not an efficient way to do this convolution
                           but a handy procedure that handles 3c data simply*/
                        if(filter_stack) 
                            result=sparse_convolve(filter_wavelet,result);
			double stalat=result.get_double("site.lat");
			double stalon=result.get_double("site.lon");
			double staelev=result.get_double("site.elev");
			double ttime=hcen.phasetime(rad(stalat),rad(stalon),
				staelev,phase);
			result.put("dir",dir);
                        dfile=dfile_base+ss.str()+".w";
                        cout << "Writing to dfile="<<dfile<<endl;
			result.put("dfile",dfile);
			/* result is in relative time.  This makes 
			result arrival time reference = predicted time  for
			this phase */
			double atime=hcen.time+ttime;
			result.rtoa(atime);
			result.put("arrival.time",atime);
			/* This program writes fake, absolute times */
			result.put("timetype","a");
			result.put("wfprocess.algorithm","RFeventstacker");
			/* Save data to both wfdisc and wfprocess */
			int irec;
			irec=dbsave(result,dbwfdisc.db,string("wfdisc"),
				mdwfdisc,am,chanmap,true);
                        // Change dfile name to split the wfdisc and wfprocess
                        // files
                        dfile=dfile+"3c";
                        result.put("dfile",dfile);
			irec=dbsave(result,dbwfprocess.db,string("wfprocess"),
				mdwfprocess,am);

			dbwfprocess.db.record=irec;
			int pwfid=dbwfprocess.get_int("pwfid");
			if(fold>1)dbstackstats.put("pwfid",pwfid);
			dbevlink.append();
			dbevlink.put("evid",evid);
			dbevlink.put("pwfid",pwfid);
			sta=result.get_string("sta");
			dbsclink.append();
			dbsclink.put("sta",sta);
			dbsclink.put("chan","3C");
			dbsclink.put("pwfid",pwfid);
			/* We also to have to save arrival and assoc rows */
			long int  arid=save_arrival(dbarrival,result,
					arrivalchan,phase);
			save_assoc(dbassoc,result,orid,arid,phase,hcen);

		}
	}
	catch (SeisppError& serr)
	{
		serr.log_error();
	}
	catch (exception& e)
	{
		cerr << "stdlib excecption caught in main.  Error message: "
			<< e.what()<<endl;
	}
	catch (...)
	{
		cerr << "Something threw an unexpected exception"<<endl;
	}
}
