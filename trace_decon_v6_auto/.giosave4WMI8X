#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <complex>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <omp.h>
#include "perf.h"
#include "DeconOperator.h"
#include "SeisppKeywords.h"
#include "TimeSeries.h"
#include "ThreeComponentSeismogram.h"
#include "ensemble.h"
#include "Metadata.h"
#include "Hypocenter.h"
#include "dbpp.h"
#include "filter++.h"
#include "resample.h"
#include "seispp.h"
#include "stack.h"

#include "SignalToNoise.h"

/*
This the deconvolution program working on conventional receiver function deconvolution,
which means taking the vertical component as the wavelet to deconvolve with the horizontal
components. Input and output are all based on antelope database. All parameters are defined
in the trace_decon.pf.
*/

/*
Modification history:
Feb, 2015 by Xiaotao Yang
	1. add functionality to calculate the trace SNR after rotation and before deconvolution.
		however, rotation is optional. SNR attribute will be added in the decon table as 
		"rawsnr" attribute.

*/

using namespace std;
using namespace SEISPP;

bool SEISPP::SEISPP_verbose(true);
/*! \brief Generic algorithm to load arrival times from a database.

In passive array processing a very common need is to extract time windows
around marked phase pick or a theoretical arrival time.  For measured
arrival times the css database has an awkward link to waveforms that 
causes problems when dealing with continuous data.  This is one element of
a set of functions aimed at dealing with this problem.  

This particular procedure aims to take an input ensemble of data and 
find matches for arrival times in an external database.  For each member
of the ensemble with a matching arrival it posts the arrival time to the
generalized header (Metadata) of the parent object.  To avoid testing 
for valid arrivals in the ensemble after this procedure is run a index
of data with valid arrivals is returned as an stl list container of ints.

\param dat is the input ensemble passed as a reference to the vector 
	container of data objects.  The type of the objects in the container
	is generic except class T MUST be an object that inherits Metadata.
	At time of this writing T could be TimeSeries, ThreeComponentSeismogram,
	or ComplexTimeSeries. 
\param dbh match handle to a Datascope database containing arrival. 
	Normally this should reference the standard catalog view.  That
	is:  event->origin->assoc->arrival subsetted to orid==prefor.
	The algorithm actually only cares that a find using the Metadata
	area of the vector data components will provide a unique match
	into the table the handle references (normally arrival).  This 
	tacitly assumes the handle has been properly constructed and the
	proper attributes have been loaded with the data to provide a unique
	match into arrival.  With css this REQUIRES that each component of
	dat MUST have had evid and/or orid posted before calling this function.  There
	is no other unambiguous way to match waveforms in this context to
	a unique arrival. 
\param keyword is the attribute name used to posted the arrival time data.

\return list of members of input vector with valid arrival times posted.
	The contents of the list can be traversed and used with operator[]
	to extract valid data from the ensemble.  
*/
list<long> LoadArrivalTimes(vector<ThreeComponentSeismogram>& dat,
                DatascopeMatchHandle& dbh,
		    const string keyword)
{
        std::vector<ThreeComponentSeismogram>::iterator d;
	long i;
	list<long> data_with_arrivals;
	const string base_error("Warning (LoadArrivalTimes): ");
	for(d=dat.begin(),i=0;d!=dat.end();++d,++i)
	{
		double atime;  //  arrival time.  
		if(d->live)
		{
		// First see if there is an arrival for this
		// station.  If not, skip it. 
			list<long> records
				=dbh.find(dynamic_cast<Metadata&>(*d),false);
			// if no arrival silently skip data for this station
			if(records.size()<=0) continue;
			if(records.size()>1)
			{
				string sta=d->get_string("sta");
				cerr << base_error 
					<< "found "
					<< records.size()
					<< " arrivals for station "
					<< sta <<endl
					<< "Using first found "
					<< "in database view"<<endl;
			}
			Dbptr db=dbh.db;
			// tricky usage here.  begin() returns
			// an iterator so the * operator gets the
			// value = record number of the match
			db.record=*(records.begin());
			char csta[10];
			if(dbgetv(db,0,"arrival.time",&atime,"sta",csta,0)
				== dbINVALID) 
			{
				string sta=d->get_string("sta");
				cerr << base_error
					<< "dbgetv failed"
					<< " in attempt to obtain"
					<< " arrival time for station"
					<< sta << endl
					<< "Data from this station"
					<< " will be dropped"<<endl;	
			}
			d->put(keyword,atime);
			data_with_arrivals.push_back(i);
		}
	}
	return(data_with_arrivals);
} 
ThreeComponentEnsemble *BuildRegularGather(ThreeComponentEnsemble& raw,
	DatascopeMatchHandle& dbh, 
	ResamplingDefinitions& rdef,
	double target_dt,
	TimeWindow processing_window)
{
	const string arrival_keyword("arrival.time");
	const double samprate_tolerance(0.01);  // fractional sample rate tolerance
	int nmembers=raw.member.size();
	auto_ptr<TimeSeries> x1,x2,x3;
	ThreeComponentEnsemble *result;
	result = new ThreeComponentEnsemble(raw);
	// An inefficiency here, but this allow us to discard dead
	// traces and problem data from ensemble as we assemble the
	// new one.
	result->member.clear();
	result->member.reserve(raw.member.size());
	// Load arrivals from database.  List returned is index into raw of
	// data with valid arrivals loaded
	list<long> data_with_arrivals;
        data_with_arrivals=LoadArrivalTimes(raw.member,dbh,arrival_keyword);
	list<long>::iterator index;
	for(index=data_with_arrivals.begin();index!=data_with_arrivals.end();++index)
	{
		ThreeComponentSeismogram d=raw.member[*index];
		if(d.live)
		{
		try {
cout << d.get_string("sta")<<" has arrival time ="
	<<strtime(d.get_double("arrival.time"))<<endl;
			d.rotate_to_standard();	
			// partial clone used to hold result
			ThreeComponentSeismogram d3c(d);  
			x1=auto_ptr<TimeSeries>(ExtractComponent(d,0));
			x2=auto_ptr<TimeSeries>(ExtractComponent(d,1));
			x3=auto_ptr<TimeSeries>(ExtractComponent(d,2));
			// resample if necessary.  Using auto_ptr to avoid temporary pointer
			// and as good practice to avoid memory leaks
			if( (abs( (d.dt)-target_dt)/target_dt) > samprate_tolerance)
			{
				*x1=ResampleTimeSeries(*x1,rdef,target_dt,false);
				*x2=ResampleTimeSeries(*x2,rdef,target_dt,false);
				*x3=ResampleTimeSeries(*x3,rdef,target_dt,false);
			}
			// This procedure returns an auto_ptr.  An inconsistency in
			// SEISPP due to evolutionary development
			x1=ArrivalTimeReference(*x1,arrival_keyword,processing_window);
			x2=ArrivalTimeReference(*x2,arrival_keyword,processing_window);
			x3=ArrivalTimeReference(*x3,arrival_keyword,processing_window);
			// safer than using attribute ns in x1
			// assumes all three components are equal length,
			// which is pretty much guaranteed here
			int ns=x1->s.size();
			d3c.ns=ns;
			d3c.dt=x1->dt;
			d3c.t0=x1->t0;
			d3c.tref=x1->tref;
			d3c.u=dmatrix(3,ns);
			// Using blas here for speed
			dcopy(ns,&(x1->s[0]),1,d3c.u.get_address(0,0),3);
			dcopy(ns,&(x2->s[0]),1,d3c.u.get_address(1,0),3);
			dcopy(ns,&(x3->s[0]),1,d3c.u.get_address(2,0),3);
			result->member.push_back(d3c);
		} catch (SeisppError serr)
		{
			// Minor maintenance issue here.  Frozen name is
			// assumed to be in Metadata area.  Avoiding second
			// handler for possible MetadataError intentionally
			string sta=d.get_string("sta");
			cerr << "Problem assembling 3C seismogram for station "
				<< sta <<endl;
			//serr.log_error();
			raw.member[*index].live=false;
			cerr << "Data for this station dropped"<<endl;
		}
		}
			
	}
	return(result);
}
void PostEvid(ThreeComponentEnsemble *d,int evid)
{
	vector<ThreeComponentSeismogram>::iterator dptr;
	for(dptr=d->member.begin();dptr!=d->member.end();++dptr)
		dptr->put("evid",evid);
}
/*! \brief Builds a standard catalog view from a CSS3.0 database.

Much passive array processing is built around the css3.0 schema.
The css3.0 schema has a standard view that defines the definitive
catalog for a network.  This is formed by the join of:
	event->origin->assoc->arrival
and is usually (as here) subsetted to only keep rows with
orid equal to the "preferred origin" (prefor).  This procedure
builds a handle to this database view.

\param dbh is a handle to the parent Datascope database 
	from which the view is to be constructed.  It need
	only reference the correct database. 
\return new handle that refers to the "standard" catalog view
	described above.
*/
DatascopeHandle StandardCatalogView(DatascopeHandle& dbh)
{
	DatascopeHandle result(dbh);
	dbh.lookup("event");
	dbh.natural_join("origin");
	string ss_to_prefor("orid==prefor");
	dbh.subset(ss_to_prefor);
	dbh.natural_join("assoc");
	dbh.natural_join("arrival");
	return(dbh);
}

void SaveResult(DatascopeHandle& dbh,
	ThreeComponentEnsemble* gather,
		AttributeMap& amo,
			MetadataList& mdlo,
				bool use_wfdisc,
					string dir,
						string dfile,
							vector<double> &all_meta,
								string dbout)
{
	const string wfdtable("wfdisc");
	const string wfptable("wfprocess");
	//prefix component code. xiaotao yang
	const string x0_name("I0");
	const string x1_name("I1");
	const string x2_name("I2");
    DatascopeHandle dbsclink(dbh);
    DatascopeHandle dbevlink(dbh);
    DatascopeHandle dbdecon(dbh);
	
	if(use_wfdisc)
	{
/*
//Commented by Xiaotao Yang. Test for outputing decon table as a database table instead of a text file.
// Modified from the similar lines for wfprocess output.
        ofstream ofs((dbout+".decon").c_str(), ios::app);
		vector<ThreeComponentSeismogram>::iterator d;
		int i;
		for(d=gather->member.begin(), i=0;d!=gather->member.end();++d,++i)
		{
			try {
				for(int j=0;j<3;j++)
				{
					TimeSeries *data=ExtractComponent(*d,j);
					switch(j)
					{
						case 0:
							data->put(string("chan"),x0_name);
						break;
						case 1:
							data->put(string("chan"),x1_name);
						break;
						case 2:
							data->put(string("chan"),x2_name);
						break;
					}
					data->put("dir",dir);
					data->put("dfile",dfile);
					int rnum=dbsave(*data,dbh.db,wfdtable,mdlo,amo);
					if(rnum!=-1)
					{
                        rnum++;
						ofs<<rnum<<"	"
						   <<data->get_string("sta")<<"	"
						   <<data->get_string("chan")<<"	"
						   <<all_meta[i*15+j*5]<<"	"
						   <<all_meta[i*15+j*5+1]<<"	"
						   <<all_meta[i*15+j*5+2]<<"	"
						   <<all_meta[i*15+j*5+3]<<"	"
						   <<all_meta[i*15+j*5+4]<<endl;
					}
					delete data;
				}
			} catch (SeisppError& serr) {
				string sta=d->get_string("sta");
				cerr << "Error saving station "<<sta<<endl;
				serr.log_error();
			}
		}
        ofs.close();
*/
		dbdecon.lookup("decon");
		vector<ThreeComponentSeismogram>::iterator d;
		int i;
		for(d=gather->member.begin(), i=0;d!=gather->member.end();++d,++i)
		{
			try {
				for(int j=0;j<3;j++)
				{
					TimeSeries *data=ExtractComponent(*d,j);
					switch(j)
					{
						case 0:
							data->put(string("chan"),x0_name);
						break;
						case 1:
							data->put(string("chan"),x1_name);
						break;
						case 2:
							data->put(string("chan"),x2_name);
						break;
					}
					data->put("dir",dir);
					data->put("dfile",dfile);
					int rnum=dbsave(*data,dbh.db,wfdtable,mdlo,amo);
					if(rnum!=-1)
					{
						rnum++;
                        dbdecon.append();
                        dbdecon.put("pwfid",rnum);
                        dbdecon.put("sta",d->get_string("sta"));
                        //cout<<"sta = "<<d->get_string("sta")<<endl;
                        switch(j)
						{
							case 0:
								dbdecon.put("chan",x0_name);
								//cout<<"rawsnr0 = "<<d->get_double("rawsnr0")<<endl;
								dbdecon.put("rawsnr",d->get_double("rawsnr0"));
							break;
							case 1:
								dbdecon.put("chan",x1_name);
								//cout<<"rawsnr1 = "<<d->get_double("rawsnr1")<<endl;
								dbdecon.put("rawsnr",d->get_double("rawsnr1"));
							break;
							case 2:
								dbdecon.put("chan",x2_name);
								//cout<<"rawsnr2 = "<<d->get_double("rawsnr2")<<endl;
								dbdecon.put("rawsnr",d->get_double("rawsnr2"));
							break;
						}
                        dbdecon.put("niteration",   int(all_meta[i*15+j*5]));
                        dbdecon.put("nspike",       int(all_meta[i*15+j*5+1]));
                        dbdecon.put("epsilon",      all_meta[i*15+j*5+2]);
                        dbdecon.put("peakamp",      all_meta[i*15+j*5+3]);
                        dbdecon.put("averamp",      all_meta[i*15+j*5+4]);
						//dbdecon.put("snr", 			1.0);
					}
					delete data;
				}
			} catch (SeisppError& serr) {
				string sta=d->get_string("sta");
				cerr << "Error saving station "<<sta<<endl;
				//serr.log_error();
			}
		}
	}
	else
	{
		dbsclink.lookup("sclink");
		dbevlink.lookup("evlink");
        dbdecon.lookup("decon");
		vector<ThreeComponentSeismogram>::iterator d;
		int i;
		for(d=gather->member.begin(),i=0;d!=gather->member.end();++d,++i)
		{
			try {
				d->put("dir",dir);
				d->put("dfile",dfile);
				d->put("timetype",string("a"));
				d->put("wfprocess.algorithm",string("trace_decon"));
				int rnum;
				if((rnum=dbsave(*d,dbh.db,wfptable,mdlo,amo))!=-1)
				{
                    rnum++;
					for(int j=0;j<3;j++)
					{
                        dbdecon.append();
                        dbdecon.put("pwfid",rnum);
                        dbdecon.put("sta",d->get_string("sta"));
                        switch(j)
						{
							case 0:
                                dbdecon.put("chan",x0_name);
                                //cout<<"rawsnr0 = "<<d->get_double("rawsnr0")<<endl;
								dbdecon.put("rawsnr",d->get_double("rawsnr0"));
                                break;
							case 1:
                                dbdecon.put("chan",x1_name);
                                //cout<<"rawsnr1 = "<<d->get_double("rawsnr1")<<endl;
								dbdecon.put("rawsnr",d->get_double("rawsnr1"));
                                break;
							case 2:
                                dbdecon.put("chan",x2_name);
                                //cout<<"rawsnr2 = "<<d->get_double("rawsnr2")<<endl;
								dbdecon.put("rawsnr",d->get_double("rawsnr2"));
                                break;
						}
                        dbdecon.put("niteration",   int(all_meta[i*15+j*5]));
                        dbdecon.put("nspike",       int(all_meta[i*15+j*5+1]));
                        dbdecon.put("epsilon",      all_meta[i*15+j*5+2]);
                        dbdecon.put("peakamp",      all_meta[i*15+j*5+3]);
                        dbdecon.put("averamp",      all_meta[i*15+j*5+4]);
					}
					dbevlink.append();
					dbevlink.put("evid",d->get_int("evid"));
					dbevlink.put("pwfid",rnum);
					dbsclink.append();
					dbsclink.put("sta",d->get_string("sta"));
					dbsclink.put("chan","3C");
					dbsclink.put("pwfid",rnum);
				}
			} catch (SeisppError& serr) {
				string sta=d->get_string("sta");
				cerr << "Error saving station "<<sta<<endl;
				//serr.log_error();
			}
		}
	}
}
void TestNaN(ThreeComponentSeismogram& dat)
{
	for(int i=0;i<dat.u.rows();i++)
		for(int j=0;j<dat.u.columns();j++)
		{
			if(isnan(dat.u(i,j)))
				{
					dat.live=false;
					cout<<"NaN found in "
						<<dat.get_string("sta")
						<<endl;
					return;
				}
		}
}
/*! \brief Multiple stage TimeInvariantFilter operator.

Sometimes it is necessary to apply a series of stages of
filtering to equalize different data sets or just to 
simply the process of defining a multiple filter chain.
This object simplifies that process by allowing the definition
of a chain of filters and a set of supplied methods to apply
these filters to data.
*/
class MultiStageFilter
{
public:
	/*! \brief Default constructors.  

	Loads a null filter definition.  That is the default is
	a do nothing operator. */
	MultiStageFilter();
	/*! \brief Construct the filter definitions from a string.

	This is currently the prime constructor.  A string is parsed
	into tokens that describe a series of filters that will be
	chained together.  The strings that define each individual
	filter type are parsed into blocks determined by the separator
	argument.  The parsed strings are currently used to construct
	a set of TimeInvariantFilter processing objects.
	An example helps explain how this would be used.  If we
	passed "DEMEAN; BW 0.5 5 2 5" and define ";" as the separator
	this would yield a two stage filter:  demean followed by a 
	0.5 to 2 Hz bandpass filter.

	\param filterlist contains the set of filter recipes to use.
	\param separator is the string used as a separator between
		the recipes for each filter description.
	*/
	MultiStageFilter(string filterlist,string separator);
	template <class T> void  apply(T& d);
private:
	list<TimeInvariantFilter> stages;
};
MultiStageFilter::MultiStageFilter()
{
	TimeInvariantFilter f(string("none"));
	stages.push_back(f);
}
MultiStageFilter::MultiStageFilter(string filterlist,string separator)
{
    try {
	const string white(" \t\n");
	int current=0,end_current;
	string stmp;
	string filterparam;
	// Strip any leading white space
	stmp=filterlist;
	if((current=stmp.find_first_not_of(white,0)) != 0)
	{
		stmp.erase(0,current);
		current=0;
	}
	int endstmp=stmp.size();
	do {
		end_current=stmp.find_first_of(separator,current);
		if(end_current<0)
		{
			filterparam.assign(stmp,current,endstmp);
			stages.push_back(TimeInvariantFilter(filterparam));
			break;
		}			
		filterparam.assign(stmp,current,end_current-current);
		stages.push_back(TimeInvariantFilter(filterparam));
		current=stmp.find_first_not_of(separator,end_current);
		current=stmp.find_first_not_of(white,current);
	}
	while(current<endstmp && current>=0);
    } catch (...) {throw;};
}
		
	
	
/* For now this only works on ensembles.  */
template <class T> void MultiStageFilter::apply(T& d)
{
    try {
	list<TimeInvariantFilter>::iterator filt;
	for(filt=stages.begin();filt!=stages.end();++filt)
		FilterEnsemble(d,*filt);
	/* Note this template could be made to work with TimeSeries
	 or ThreeComponentSeismogram components individually if 
	we used a typeid check on T and then used this line
	for that type of beast:  filt->apply(d);
	*/
    } catch (...) {throw;};
}
void ApplyFST(ThreeComponentSeismogram& e,Hypocenter& hypo,string type)
{
	double vp0(6.0),vs0(3.5);
	double lat,lon,elev;
	lat=e.get_double("lat");
	lon=e.get_double("lon");
	elev=e.get_double("elev");
	SlownessVector u=hypo.pslow(lat,lon,elev);
	if(type=="zrt")
	{
		// First the horizonal rotation
		SphericalCoordinate scor;
		scor.phi=atan2(u.uy,u.ux);
		scor.theta=0.0;
		scor.radius=1.0;
		// after this transformation x1=transverse horizontal
		// x2=radial horizonal, and x3 is still vertical
		e.rotate(scor);
	}
	else if(type=="fst")
		e.free_surface_transformation(u,vp0,vs0);
	else
		cout<<"No rotation applied before deconvolution!"<<endl;
	
}
void ApplyFST(ThreeComponentEnsemble& e,Hypocenter& hypo,string type)
{
	for(int i=0;i<e.member.size();++i)
		ApplyFST(e.member[i],hypo,type);
	/*double vp0(6.0),vs0(3.5);
	if(type=="zrt")
		for(int i=0;i<e.member.size();++i)
		{
			double lat,lon,elev;
			lat=e.member[i].get_double("lat");
			lon=e.member[i].get_double("lon");
			elev=e.member[i].get_double("elev");
			SlownessVector u=hypo.pslow(lat,lon,elev);
		    // First the horizonal rotation
		    SphericalCoordinate scor;
		    scor.phi=atan2(u.uy,u.ux);
		    scor.theta=0.0;
		    scor.radius=1.0;
		    // after this transformation x1=transverse horizontal
		    // x2=radial horizonal, and x3 is still vertical
		    e.member[i].rotate(scor);
		}
	else if(type=="fst")
		for(int i=0;i<e.member.size();++i)
		{
			double lat,lon,elev;
			lat=e.member[i].get_double("lat");
			lon=e.member[i].get_double("lon");
			elev=e.member[i].get_double("elev");
			SlownessVector u=hypo.pslow(lat,lon,elev);
			e.member[i].free_surface_transformation(u,vp0,vs0);
		}*/
}
void ApplyKills(ThreeComponentEnsemble *d, set<int>& kills)
{
    set<int>::iterator kptr;
    int nmembers=d->member.size();
    for(kptr=kills.begin();kptr!=kills.end();++kptr)
    {
        int i=(*kptr);
        if((i<0) || (i>=nmembers) ) throw i;
        d->member[i].live=false;
    }
}
unsigned int nextPowerOf2(unsigned int n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}
void usage()
{
	cerr << "Version: 2015.12.01" <<endl;
	cerr << "trace_decon db [-d output_dir -o outputdb -n start_evid -pf pfname]" << endl;
	exit(-1);
}

int main(int argc, char **argv)
{
	// This is a C++ thing.  Not required on all platforms, but
	// recommended by all the books.  
	ios::sync_with_stdio();
	// This is a standard C construct to crack the command line
	if(argc<2) usage();
	string pfin(argv[0]);
	string dbin(argv[1]);
	string dbout=dbin;
	string output_dir="tracedecon";
	string stnumstring("NULL");
	for(int i=2;i<argc;++i)
	{
		if(!strcmp(argv[i],"-V"))
			usage();
		else if(!strcmp(argv[i],"-pf"))
		{
			++i;
			if(i>=argc) usage();
			pfin = string(argv[i]);
		}
		else if(!strcmp(argv[i],"-o"))
		{
			++i;
			if(i>=argc) usage();
			dbout=string(argv[i]);
		}
		else if(!strcmp(argv[i],"-d"))
		{
			++i;
			if(i>=argc) usage();
			output_dir=string(argv[i]);
		}
		else if(!strcmp(argv[i],"-n"))
		{
			++i;
			if(i>=argc) usage();
			stnumstring=string(argv[i]);
		}
		else
			usage();
	}
	// Standard Antelope C construct to read a parameter file
	// Well almost standard.  const_cast is a C++ idiom required
	// due to a collision of constantness with Antelope libraries.
	Pf *pf;
        if(pfread(const_cast<char *>(pfin.c_str()),&pf)) 
	{
		cerr << "pfread error for pf file="<<pfin<<".pf"<<endl;
		exit(-1);
	}
	Metadata control(pf);
	MetadataList mdens=pfget_mdlist(pf,"Ensemble_mdlist");
	MetadataList mdtrace=pfget_mdlist(pf,"station_mdlist");
	MetadataList mdlo;
	try {
		// Get set of control parameters
		string netname=control.get_string("netname");
		string phase=control.get_string("phase");
		double ts,te;
		ts=control.get_double("data_window_start");
		te=control.get_double("data_window_end");
		TimeWindow datatwin(ts,te);
		double tpad=control.get_double("data_time_pad");
		ts=control.get_double("processing_window_start");
		te=control.get_double("processing_window_end");
		TimeWindow processing_twin(ts,te);
		ts=control.get_double("noise_window_start");
		te=control.get_double("noise_window_end");
		TimeWindow noise_twin(ts,te);
		StationChannelMap stachanmap(pf);
		string schemain=control.get_string("InputAttributeMap");
		string schemaout=control.get_string("OutputAttributeMap");
		double target_dt=control.get_double("target_sample_interval");
		int time_shift=-control.get_double("processing_window_start")/target_dt;
		if(time_shift<0) time_shift=0;
		control.put("shaping_wavelet_dt",target_dt);
		control.put("sample_shift",time_shift);
		ResamplingDefinitions rdef(pf);
		
		//string output_dir=control.get_string("output_dir");
		bool use_wfdisc=control.get_bool("use_wfdisc");
		
		if(use_wfdisc) mdlo=pfget_mdlist(pf,"output_mdlist_wfdisc");
		else
			mdlo=pfget_mdlist(pf,"output_mdlist_wfprocess");
		//read in snr time window parameters. xiaotao yang
		double snr_signal_window_start=control.get_double("snr_signal_window_start");
		double snr_signal_window_end=control.get_double("snr_signal_window_end");
		double snr_noise_window_start=control.get_double("snr_noise_window_start");
		double snr_noise_window_end=control.get_double("snr_noise_window_end");
		if(fabs(snr_signal_window_end-snr_signal_window_start)<=1.0e-10 ||
			fabs(snr_noise_window_end-snr_noise_window_start)<=1.0e-10)
		{
			cerr<<"ERROR: Wrong time window parameters for snr calculation!"<<endl;
			exit(-1);
		}
		// First get all the database components assembled.
		// input data, output data
		DatascopeHandle dbh(dbin,false);
		/* Build a match handle for arrivals using the standard catalog
		view of the join of event->origin->assoc->arrival */
		AttributeMap am(schemain);  
		AttributeMap amo(schemaout);  
		DatascopeHandle dbcatalog=StandardCatalogView(dbh);
		list<string> matchkeys;
		matchkeys.push_back(string("sta"));
		matchkeys.push_back(string("evid"));
		DatascopeMatchHandle dbhm(dbcatalog,string(""),matchkeys,am);
        
		/* Build a event view to drive the process*/
		DatascopeHandle dbhv(dbh);
		dbhv.lookup("event");
		dbhv.natural_join("origin");
		dbhv.subset(string("orid==prefor"));
		
		/* Define output database*/
		DatascopeHandle dbho(dbh);
		if(dbout!=dbin)
			dbho=DatascopeHandle(dbout,false);
		dbho.lookup(string("wfprocess"));
		
		/* Need an editor window for data plotting.   Deconvolved
		   data will be replace original data. */
		Metadata dpdef(pf,string("data_plot_properties"));
		///TraceEditPlot dataplot(dpdef);
		/* these raw pointers are used to hold dynamic ensembles
		   below.  Usual warning about pointers. */
		ThreeComponentEnsemble *rawdata=NULL, 
				*regular_gather=NULL,
				*decondata=NULL;
		SeismicArray *stations=NULL;
		/* Assorted variables needed in the loop below */
		double lat,lon,depth,otime;
		int evid;
		int record;
		string filter_param;
		SimpleDecon *decon_op;
		dbhv.rewind();
		if(stnumstring.compare("NULL")!=0)
		{
			long idst=atol(stnumstring.c_str());
			for(int i=0;i<dbhv.number_tuples();++i,++dbhv)
			{
				if(idst==dbhv.get_long(string("evid")))
					break;
			}
			++dbhv;
			cout<<"First evid to be process is "<<dbhv.get_long(string("evid"))<<endl;
		}
		long recordnum=dbhv.number_tuples()-dbhv.current_record();
        //for now, test only!!
        //for(int i=0;i<514;i++) ++dbhv;
		for(record=0;record<recordnum;++record,++dbhv)
		{
			int num_gather;
			lat=dbhv.get_double(string("origin.lat"));
			lon=dbhv.get_double(string("origin.lon"));
			depth=dbhv.get_double(string("origin.depth"));
			otime=dbhv.get_double(string("origin.time"));
			evid=dbhv.get_int(string("evid"));
			// origin coordinates are degrees in the db,
			// but must be radians for hypocenter object
			lat=rad(lat);
			lon=rad(lon);
			
			Hypocenter hypo(lat,lon,depth,otime,
				string("tttaup"),string("iasp91"));
			// Define the filter
			try {
				filter_param=control.get_string("filter");
			} catch (MetadataGetError mdge) {
				cerr	<<"filter attribute not defined"
					<<" using default of DEMEAN"<<endl;
				filter_param=string("DEMEAN");
			}
			MultiStageFilter* filt=new MultiStageFilter(filter_param,string(";"));
			// On the first record we need to load the station
			// geometry object
			if(record==0)
			{
				stations = new SeismicArray(
					dynamic_cast<DatabaseHandle&>(dbh),
					hypo.time,netname);
			}
			else
			{
				TimeWindow test(hypo.time,hypo.time+2000.0);
				if(!stations->GeometryIsValid(test))
				{
					delete stations;
					stations = new SeismicArray(
					dynamic_cast<DatabaseHandle&>(dbh),
					hypo.time,netname);
				}
			}
			// Read the raw data using the time window based constructor
			try{
				rawdata=array_get_data(*stations,hypo,
					phase,datatwin,tpad,dynamic_cast<DatabaseHandle&>(dbh),
					stachanmap,mdens,mdtrace,am);
			}catch (SeisppError& serr)
			{
				cout<<"bad event"<<endl;
				//serr.log_error();
				continue;
			}catch (...)
			{
				cout<<"Unexpected error from array_get_data!"<<endl;
				continue;
			}
			if(rawdata->member.size()==0)
			{
			cout<<"bad event with problem wf data"<<endl;
			continue;
			}
			ThreeComponentEnsemble *rawdata_filt;
			try{
				///do
				///{
					rawdata_filt=new ThreeComponentEnsemble(*rawdata);
					// Filter the raw data.
					filt->apply(*rawdata_filt);
					PostEvid(rawdata_filt,evid);

					num_gather=rawdata_filt->member.size();
					for(int i=0;i<num_gather;i++)
					{
						double lat=stations->array[rawdata_filt->member[i].get_string("sta")].lat;
						double lon=stations->array[rawdata_filt->member[i].get_string("sta")].lon;
						double elev=stations->array[rawdata_filt->member[i].get_string("sta")].elev;
						rawdata_filt->member[i].put("lat",lat);
						rawdata_filt->member[i].put("lon",lon);
						rawdata_filt->member[i].put("elev",elev);
						//cout<<rawdata->member[i].get_string("sta")<<":"<<lat<<":"<<lon<<":"<<elev<<endl;
					}

					regular_gather=BuildRegularGather(*rawdata_filt, dbhm,rdef,target_dt,
							processing_twin);
					for(int i=0;i<num_gather;i++)
					{
						if(regular_gather->member[i].ns<time_shift)
						{
							cout<<"bad event with problem waveform that's too short"<<endl;
							throw i;
						}
					}
						    
					if(control.get_bool("apply_rotation"))
					{
						string rotation_type=control.get_string("rotation_type");
						ApplyFST(*regular_gather,hypo,rotation_type);
						//cout<<"applied FST"<<endl;
					}
						// regular gather is now assumed to contain data with
						// a common start time.  
						//regular_gather->put(string("moveout"),0.0);
	
						//cout<<"Number of Samples:"<<regular_gather->member[0].ns<<endl;
						//cout<<nextPowerOf2(regular_gather->member[0].ns)<<endl;
						///dataplot.plot(*regular_gather,true);
		
					string filt_str;
						///cout<<"change filter parameters?"<<endl;
			    ///std::getline(cin,filt_str);
						///if(filt_str=="no")
						///	break;
						///else
						///{
						///	delete filt;
						///	filt=new MultiStageFilter(filt_str,string(";"));
			    ///    delete rawdata_filt;
						///}
				///}while(1);
			}catch (...)
			{
				cout<<"bad event with problem waveform"<<endl;
				//serr.log_error();
				continue;
			}
			delete filt;
			delete rawdata;
			rawdata=new ThreeComponentEnsemble(*rawdata_filt);
			delete rawdata_filt;
            
			num_gather=regular_gather->member.size();
			int maxns=0;
			for(int i=0;i<num_gather;i++)
			{
				if(regular_gather->member[i].ns>maxns)
					maxns=regular_gather->member[i].ns;
			}
			if(maxns<=0)
			{
				cout<<"bad event that has a number of sample eq or lt 0"<<endl;
				continue;
			}
			
			control.put("operator_nfft",static_cast<int>(nextPowerOf2(maxns)));
			control.put("taper_length",maxns);
			string decon_type=control.get_string("deconvolution_type");
			
			if(decon_type=="least_square")
				decon_op=new SimpleLeastSquareDecon(control);
			else if(decon_type=="water_level")
				decon_op=new SimpleWaterLevelDecon(control);
			else if(decon_type=="multi_taper")
				decon_op=new SimpleMultiTaperDecon(control);
			else if(decon_type=="iterative")
				decon_op=new SimpleGeneralIterDecon(control);
			else
			{
				cout<<"WARNING: deconvolution_type is unavailable, set to least_square by default."<<endl;
				decon_type="least_square";
				decon_op=new SimpleLeastSquareDecon(control);
			}
			
			
			decondata=new ThreeComponentEnsemble(*regular_gather);
			
			//for debug
			//omp_set_num_threads(1);
			
			vector<vector<double> > all_data;
			all_data.reserve(num_gather*3);
			
			for(int i=0;i<num_gather;i++)
			{
				for(int j=0;j<3;j++)
				{
					TimeSeries *data=ExtractComponent(regular_gather->member.at(i),j);
					all_data.push_back(data->s);
					delete data;
				}
			}

			bool multi_taper_on=false;
			if(decon_type=="multi_taper" ||
				(decon_type=="iterative" && control.get_string("iterative_decon_type")=="multi_taper"))
				multi_taper_on=true;
            
			vector<vector<double> > all_noise;
			if(multi_taper_on)
			{
				const string arrival_keyword("arrival.time");
				all_noise.reserve(num_gather);
				list<long> data_with_arrivals;
				data_with_arrivals=LoadArrivalTimes(rawdata->member,dbhm,arrival_keyword);
				list<long>::iterator index;
				for(index=data_with_arrivals.begin();index!=data_with_arrivals.end();++index)
				{
					ThreeComponentSeismogram d=rawdata->member[*index];
					if(d.live)
					{
						d.rotate_to_standard();
						if(control.get_bool("apply_rotation"))
						{
							string rotation_type=control.get_string("rotation_type");
							ApplyFST(d,hypo,rotation_type);
						}
						ThreeComponentSeismogram d3c(d);
						auto_ptr<TimeSeries> noise;
						noise=auto_ptr<TimeSeries>(ExtractComponent(d,2));
						if( (abs( (noise->dt)-target_dt)/target_dt) > 0.01)
						{
							*noise=ResampleTimeSeries(*noise,rdef,target_dt,false);
						}
						noise=ArrivalTimeReference(*noise,arrival_keyword,noise_twin);
						all_noise.push_back(noise->s);
						noise.release();
					}
				}
			}
			
			vector<double> all_meta;
			if(decon_type=="iterative")
			{
				all_meta.reserve(num_gather*3*5);
				all_meta.resize(num_gather*3*5);
				
				// ADD SNR ATTRIBUTE IN DECON TABLE.
				//all_meta.reserve(num_gather*3*6);
				//all_meta.resize(num_gather*3*6);
			}
			
			vector<vector<double> > all_result;
			all_result.reserve(num_gather*3);
			all_result.resize(num_gather*3);
			#pragma omp parallel
			{
			#pragma omp for schedule(guided, 1)
				for(int i=0;i<num_gather;i++)
				{
					SimpleDecon *para_decon_op;
					if(decon_type=="least_square")
						para_decon_op=new SimpleLeastSquareDecon(*dynamic_cast<SimpleLeastSquareDecon *>(decon_op));
					else if(decon_type=="water_level")
						para_decon_op=new SimpleWaterLevelDecon(*dynamic_cast<SimpleWaterLevelDecon *>(decon_op));
					else if(decon_type=="multi_taper")
						para_decon_op=new SimpleMultiTaperDecon(*dynamic_cast<SimpleMultiTaperDecon *>(decon_op));
					else if(decon_type=="iterative")
						para_decon_op=new SimpleGeneralIterDecon(*dynamic_cast<SimpleGeneralIterDecon *>(decon_op));
					else
						para_decon_op=new SimpleLeastSquareDecon(*dynamic_cast<SimpleLeastSquareDecon *>(decon_op));
					
					if(multi_taper_on)
					{
						if(decon_type=="multi_taper")
							dynamic_cast<SimpleMultiTaperDecon *>(para_decon_op)->loadnoise(all_noise[i]);
						else
							dynamic_cast<SimpleGeneralIterDecon *>(para_decon_op)->loadnoise(all_noise[i]);
					}
						
					para_decon_op->load(all_data[i*3+2],all_data[i*3]);
					all_result[i*3]=para_decon_op->getresult();
					cout<<"THREAD:"<<omp_get_thread_num()<<" processed "
						<<regular_gather->member.at(i).get_string("sta")
						<<":0"<<endl;
					
					if(decon_type=="iterative")
					{
						int spikecount=0;
						double peakamp=0;
						double rms=0;
						for(int j=0;j<all_result[i*3].size();j++)
						{
							if(abs(all_result[i*3][j])>1e-15)
								spikecount++;
							if(abs(all_result[i*3][j])>peakamp)
								peakamp=abs(all_result[i*3][j]);
							rms+=all_result[i*3][j]*all_result[i*3][j];
						}
						all_meta[i*3*5]=dynamic_cast<SimpleGeneralIterDecon *>(para_decon_op)->numberiter()+1;
						all_meta[i*3*5+1]=spikecount;
						all_meta[i*3*5+2]=dynamic_cast<SimpleGeneralIterDecon *>(para_decon_op)->epsilon();
						all_meta[i*3*5+3]=peakamp;
						all_meta[i*3*5+4]=sqrt(rms/all_result[i*3].size());
					}

					para_decon_op->load(all_data[i*3+2],all_data[i*3+1]);
					all_result[i*3+1]=para_decon_op->getresult();
					cout<<"THREAD:"<<omp_get_thread_num()<<" processed "
						<<regular_gather->member.at(i).get_string("sta")
						<<":1"<<endl;
					
					if(decon_type=="iterative")
					{
						int spikecount=0;
						double peakamp=0;
						double rms=0;
						for(int j=0;j<all_result[i*3+1].size();j++)
						{
							if(abs(all_result[i*3+1][j])>1e-15)
								spikecount++;
							if(abs(all_result[i*3+1][j])>peakamp)
								peakamp=abs(all_result[i*3+1][j]);
							rms+=all_result[i*3+1][j]*all_result[i*3+1][j];
						}
						all_meta[(i*3+1)*5]=dynamic_cast<SimpleGeneralIterDecon *>(para_decon_op)->numberiter()+1;
						all_meta[(i*3+1)*5+1]=spikecount;
						all_meta[(i*3+1)*5+2]=dynamic_cast<SimpleGeneralIterDecon *>(para_decon_op)->epsilon();
						all_meta[(i*3+1)*5+3]=peakamp;
						all_meta[(i*3+1)*5+4]=sqrt(rms/all_result[i*3+1].size());
					}
					
					para_decon_op->load(all_data[i*3+2],all_data[i*3+2]);
					all_result[i*3+2]=para_decon_op->getresult();
					cout<<"THREAD:"<<omp_get_thread_num()<<" processed "
						<<regular_gather->member.at(i).get_string("sta")
						<<":2"<<endl;
					
					if(decon_type=="iterative")
					{
						int spikecount=0;
						double peakamp=0;
						double rms=0;
						for(int j=0;j<all_result[i*3+2].size();j++)
						{
							if(abs(all_result[i*3+2][j])>1e-15)
								spikecount++;
							if(abs(all_result[i*3+2][j])>peakamp)
								peakamp=abs(all_result[i*3+2][j]);
							rms+=all_result[i*3+2][j]*all_result[i*3+2][j];
						}
						all_meta[(i*3+2)*5]=dynamic_cast<SimpleGeneralIterDecon *>(para_decon_op)->numberiter()+1;
						all_meta[(i*3+2)*5+1]=spikecount;
						all_meta[(i*3+2)*5+2]=dynamic_cast<SimpleGeneralIterDecon *>(para_decon_op)->epsilon();
						all_meta[(i*3+2)*5+3]=peakamp;
						all_meta[(i*3+2)*5+4]=sqrt(rms/all_result[i*3+2].size());
					}
					
					delete para_decon_op;
				}
			}
			
			delete rawdata;
			
			//define timewindows for signal and noise ratio calculation.
			//Xiaotao Yang
			TimeWindow noise(snr_noise_window_start,snr_noise_window_end),
						signal(snr_signal_window_start,snr_signal_window_end);
			// SNR_rms for ThreeComponentSeismogram is not ready yet. Use TimeSeries version for now.
			//ensemble_SNR_rms<ThreeComponentEnsemble, ThreeComponentSeismogram>(*regular_gather,signal,noise,SEISPP::snr_keyword);

			for(int i=0;i<num_gather;i++)
			{
				vector<TimeSeries> decon_result_vector;
				for(int j=0;j<3;j++)
				{
					double snr_rms;
					TimeSeries *data=ExtractComponent(regular_gather->member.at(i),j);
					//get SNR before replacing the raw data with the deconvolved data. XIAOTAO YANG
					//cout<<"t0 in timeseries = "<<data->t0<<endl;
					snr_rms=SNR_rms(*data,signal,noise);
					//data->put(SEISPP::snr_keyword,snr_rms);

					data->s=all_result[i*3+j];
					decon_result_vector.push_back(*data);
					//CHECK SNR VALUES. xiaotao yang
					//DEBUG
					/*
					cout<<"sta = "<<data->get_string("sta")<<", comp = "
						<<j<<", snr = "
						<<decon_result_vector[j].get_double(SEISPP::snr_keyword)<<endl;
					*/
					//put snr to decondata. Xiaotao Yang
					switch(j)
					{
						case 0:
							decondata->member[i].put("rawsnr0",snr_rms);
							break;
						case 1:
							decondata->member[i].put("rawsnr1",snr_rms);
							break;
						case 2:
							decondata->member[i].put("rawsnr2",snr_rms);
							break;
					}
					
					delete data;
				}
				ThreeComponentSeismogram decon_result_seismogram(decon_result_vector,2);
				decondata->member[i].u=decon_result_seismogram.u;
				TestNaN(decondata->member[i]);
				if(control.get_bool("apply_rotation")&&control.get_bool("rotate_back"))
					decondata->member[i].rotate_to_standard();
			}

			/*cout << "Edit deconvolved data to remove bad traces"
			  <<endl
			  <<"Type X in the data window or use the menu to continue"
			  <<endl;
			*/
			///dataplot.plot(*decondata,true);
			///set<int> kill_list=dataplot.report_kills();
			/*try {
				ApplyKills(decondata,kill_list);
			}catch(int ierr)
			{
				cerr << "Fatal problem with TraceEditPlot return"
					<<endl
					<< "Editor returned invalid kill index of "
					<< ierr<<endl
					<< "Must be in range 0 to "
					<< decondata->member.size()-1<<endl;
				exit(-1);
			}*/
			//convert time reference to absolute time
			const string arrival_keyword("arrival.time");
			for(int i=0;i<num_gather;i++)
			{
				decondata->member.at(i).rtoa(regular_gather->member.at(i).get_double(arrival_keyword));
			}
			
			stringstream ss;
			ss << evid;
			string dfile = ss.str();
			SaveResult(dbho,decondata,amo,mdlo,use_wfdisc,output_dir,decon_type+"_"+dfile,all_meta,dbout);
			delete regular_gather;
			delete decondata;
			delete decon_op;
		}
	}
	catch (SeisppError& serr)
	{
		serr.log_error();
	}
	catch (...)
	{
		cerr << "Something threw an unexpected exception"<<endl;
	}
}

