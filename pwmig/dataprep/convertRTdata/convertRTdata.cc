/* Crude program to beat SAC radial receiver function into input format
for pwmig.  That means 3c data objects.  See design file in this directory
for more info. */
#include <math.h>
#include <sstream>
#include "perf.h"
#include "coords.h"
#include "stock.h"
#include "seispp.h"
#include "dbpp.h"
#include "ThreeComponentSeismogram.h"
#include "TimeSeries.h"
#include "Hypocenter.h"
#include "ArrivalUpdater.h"
using namespace std;
using namespace SEISPP;
void usage()
{
	cerr << "convertRTdata dbin dbout [-pf pffile -V]"<<endl;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	if(argc<3) usage();
	string dbin(argv[1]);
	string dbout(argv[2]);
	const string outtable("wfprocess");
	cout << "convertRTdata:  read from db="<<dbin<<" writing results to db="<<dbout<<endl;
	string pffile("convertRTdata");
	int i;
	for(i=3;i<argc;++i)
	{
		string argstr(argv[i]);
		if(argstr=="-pf") 
		{
			++i;
			if(i>=argc) usage();
			pffile=string(argv[i]);
		}
		else if(argstr=="-V")
			SEISPP_verbose=true;
		else
			usage();
	}
	Pf *pf;
	if(pfread(const_cast<char *>(pffile.c_str()),&pf))
	{
		cerr << "pfread error"<<endl;
		usage();
	}
	MetadataList mdl=pfget_mdlist(pf,"trace_mdl");
	MetadataList mdlout=pfget_mdlist(pf,"output_mdl");
	MetadataList mdlassoc=pfget_mdlist(pf,"save_assoc_mdl");
	MetadataList mdlarrival=pfget_mdlist(pf,"save_arrival_mdl");
	MetadataList mdlwfdisc=pfget_mdlist(pf,"save_wfdisc_mdl");
	vector<string>chanmap;
	chanmap.push_back("E");
	chanmap.push_back("N");
	chanmap.push_back("Z");
	Metadata control(pf);
	bool start_time_is_origin=control.get_bool("start_time_is_origin");
	/* Current sources of deconvolved data I've dealt with either have
	radial only or r and t */
	string Rchan=control.get_string("radial_channel_code");
	string Tchan=control.get_string("tangential_channel_code");
	bool save_wfdisc=control.get_bool("save_wfdisc");
        /* QC problem with EARS data found a problem with high 
           amplitude traces that dominated the stack. This parameter
           will autokill (with a message) any trace with an rms amplitude
           exceeding this size on either radial or transverse */
        double rms_kill_threshold=control.get_double("rms_kill_threshold");
		
	try {
		AttributeMap am("css3.0");
		DatascopeHandle dbhi(dbin,true);
		dbhi.lookup("event");
		dbhi.natural_join("origin");
		dbhi.natural_join("assoc");
		dbhi.natural_join("arrival");
		string sstr("iphase=~/P/");
		dbhi.subset(sstr);
		list<string> jk1,jk2;
		jk1.push_back(string("arrival.time"));
		jk1.push_back(string("sta"));
		jk2.push_back(string("wfdisc.time::wfdisc.endtime"));
		jk2.push_back(string("sta"));
		//dbhi.join("wfdisc",jk1,jk2);
		//dbhi.natural_join("wfdisc");
		dbhi.leftjoin("wfdisc",jk2,jk1);
		dbhi.natural_join("site");
		list<string> sortkeys;
		sortkeys.push_back("origin.time");
		sortkeys.push_back("sta");
		dbhi.sort(sortkeys);

		cout << "Working input view has "<<dbhi.number_tuples()
			<<" rows"<<endl;
		if(dbhi.number_tuples()<=0)
		{
			cerr << "No data to process.  Exiting.  "<<endl
				<< "Working view "
				<<"(wfdisc<-event->origin->assoc->arrival->site)"
				<< " is empty"<<endl;
			exit(-1);
		}
		DatascopeHandle dbho(dbout,false);
		dbho.lookup("wfprocess");
		/* Need to put stuff in evlink and sclink for pwmig */
		DatascopeHandle dbhosc(dbho);
		dbhosc.lookup("sclink");
		DatascopeHandle dbhoev(dbho);
		dbhoev.lookup("evlink");
		DatascopeHandle dbhwfdisc(dbho);
		dbhwfdisc.lookup("wfdisc");
		/* We use this as a convenient way to create a new assoc and arrival */
		ArrivalUpdater AUhandle(dynamic_cast<DatabaseHandle&>(dbho),mdlassoc,mdlarrival,string("css3.0"));
		int j,k;
		double slat,slon,sdepth,stime;
		string method("tttaup"),model("iasp91");
		double stalat,stalon,staelev;
		double seaz,az;
		double a,b;
		double atime;
		int evid,lastevid,pwfid;
		string sta,laststa;
		/* frozen for now */
		string outdir("wf3c");
		string dfilebase("RFdata");
		dbhi.rewind();
		laststa=dbhi.get_string("sta");
		lastevid=dbhi.get_int("evid");
		/* This holds radial and transverse - hence 2*/
		vector<TimeSeries> d;
		int ic(0);  // always 0 or 1 */
		int tracecount(0);
		int ntuples=dbhi.number_tuples();
		for(i=0;i<ntuples;++i,++dbhi)
		{
			/* This holds current time series data read.  Accumulate until sta or event change */
			TimeSeries d_read(dynamic_cast<DatabaseHandle&>(dbhi),mdl,am);
			/* ic must index member of d just pushed. Loop assumes ic=0 on first pass
			and is incremented thereafter.  Reset to 0 when new event-sta encounterd */
			sta=d_read.get_string("sta");
			evid=d_read.get_int("evid");
                        /* We compute rms here assuming we have no 
                           gaps. There is a fancier way to do this in
                           the seispp library to deal with gaps but that
                           baggage is not needed here */
                        double rms=dnrm2(d_read.s.size(),&(d_read.s[0]),1);
			if(SEISPP_verbose) cout << "Indices i,sta,chan,evid:  "
				 <<i<<", "<<sta
                                <<", "<<d_read.get_string("chan")
                                <<", "<<evid
                                <<".  RMS amplitude="<<rms<<endl;
                        if(rms>rms_kill_threshold) 
                        {
                            cout << "High rms for "<<sta<<" on channel "
                                <<d_read.get_string("chan")<<endl
                                << "Marking for deletion"<<endl;
                            d_read.live=false;
                        }
			if(sta==laststa && evid==lastevid && (i != (ntuples-1)))
			{
				if(tracecount>2)
				{
					cerr << "Warning:  irregular data  for station = "<<sta
						<<" and evid="<<evid<<endl
						<< "tracecount="<<tracecount<<" indicates duplicates"<<endl
						<< "Last in will overwrite earlier (raw db order) entries"<<endl;
				}
				d.push_back(d_read);
				++ic;
				++tracecount;
			}
			else
			{
				/* This is do so we consistently use d[0] as pattern for 3c data */
				ic=0;
				/* We need to compute an equivalent transformation matrix for
				these data from the hypocenter information.  We'll let this exit
				with an exception if origin information is not present. */
				slat=d[ic].get_double("source_lat"); // ic always 1 in this section
				slon=d[ic].get_double("source_lon");
				evid=d[ic].get_int("evid");
				cout << "Processing event "<<evid
					<< " located at (lat,lon)=("<<slat
					<<", "<<slon<<") "
					<<" for station="<<laststa<<endl;
	
				slat=rad(slat);
				slon=rad(slon);
				sdepth=d[ic].get_double("source_depth");
				stime=d[ic].get_double("source_time");
				Hypocenter h(slat,slon,sdepth,stime,method,model);
				stalat=d[ic].get_double("sta_lat");
				stalon=d[ic].get_double("sta_lon");
				staelev=d[ic].get_double("sta_elev");
				stalat=rad(stalat);
				stalon=rad(stalon);
				if(start_time_is_origin)
				{
				/* We have to fudge the start time because
				Fenglin's rf code seems to not keep absolute
				time.  Well, it might be in sac2db, but in
				any case the arrival.time of P is the same
				as the origin time.  Thus, to fix this start 
				time we find the time difference between atime
				and start time of the trace, compute P travel time,
				and then correct t0 to make P time match theoretical
				time exactly. */
				atime=d[ic].get_double("atime");
				double ptime=h.ptime(stalat,stalon,staelev);	
				double oldt0=d[ic].t0;
				double pdt=atime-oldt0;
				cout << laststa << " original atime - t0 ="<<pdt<<endl;
				cout << "origin time and atime difference="<< stime-atime<<endl;
				d[ic].t0=stime+ptime-pdt;
				cout << "t0 of this seismogram changed by "
					<< d[ic].t0-oldt0 
					<< " seconds"<<endl;
				/* correct the arrival time and we'll post it below to the new 3c trace */
				atime=stime + ptime;
				}
				else
				{
					atime=d[ic].get_double("atime");
				}
				
				seaz=h.seaz(stalat,stalon);
				/* radial is defined as propagation direction, but seaz is back azimuth */
				az=seaz-M_PI;
				if(az<(-2.0*M_PI)) az+=(2.0*M_PI);
				if(SEISPP_verbose) cout <<" station="<<laststa
					<<" computed radial azimuth="<<deg(az)
					<<" Distance(deg)="<<deg(h.distance(stalat,stalon))<<endl;
				/* Now we have to set the transformation matrix from az.  Convention here
				is we put R into x2, so the rotation angle turns out to be -az when you 
				work through the difference between rotation from the x1 axis and an
				azimuth.  In seispp teh transformation is the forward transformation so
				the signs for a and b are as list */
				a=cos(-az);
				b=sin(-az);
				/* by loading these in the metadata for d they will be cloned
				and used to construct the transformation matrix for d3c. */
				d[ic].put("U11",a);
				d[ic].put("U21",-b);
				d[ic].put("U31",0.0);
				d[ic].put("U12",b);
				d[ic].put("U22",a);
				d[ic].put("U32",0.0);
				d[ic].put("U13",0.0);
				d[ic].put("U23",0.0);
				d[ic].put("U33",1.0);
				/* this uses the Metadata from d to create a template or 3c data.*/
				ThreeComponentSeismogram d3c(dynamic_cast<Metadata&>(d[ic]),false);
				d3c.t0=d[ic].t0;
				d3c.ns=d[ic].ns;
				d3c.dt=d[ic].dt;
				d3c.components_are_orthogonal = true;
				/* initialize */
				for(k=0;k<d[ic].ns;++k)
					for(j=0;j<3;++j) d3c.u(j,k)=0.0;
				/* stuff r into 1 (x2) and t into 0 (x1).  Leave x3 0 */
				int iic,jc;
                                bool testlive;
				for(iic=0,testlive=true;iic<tracecount;++iic)
				{
                                    if(d[iic].live)
                                    {
					string chan;
					chan=d[iic].get_string("chan");
                                        if(chan.find(Rchan)!=string::npos)
						jc=1;
					else if(chan.find(Tchan)!=string::npos)
						jc=0;
					else
					{
						cerr << "Fatal: do not know how to handle channel "
							<< chan <<endl;
						exit(-1);
					}
					for(k=0;k<d[iic].ns;++k) d3c.u(jc,k)=d[iic].s[k];
                                    }
                                    else
                                    {
                                        testlive=false;
                                        break;
                                    }
				}
				if(testlive) 
                                {
                                    d3c.live=true;
        				d3c.rotate_to_standard();
        				/* We need to post these things or we're screwed */
        				d3c.put("dir",outdir);
        				char buf[128];
        				stringstream ss(buf);
        				ss<<dfilebase<<lastevid<<'\0';
        				d3c.put("dfile",ss.str());
        				d3c.put("wfprocess.algorithm","convertRF");
        				d3c.put("timetype","a");
        //cout << "3c data for sta "<<d3c.get_string("sta")<<" t0="<<strtime(d3c.t0)<<endl;
        				int rec=dbsave(d3c,dbho.db,outtable,mdlout,am);
        				if(rec<0)
        				{
					    cerr << "dbsave failed working on input record="<<i<<endl;
					    exit(-1);
				        }
				        if(save_wfdisc)
				        {
					   dbsave(d3c,dbhwfdisc.db,
						string("wfdisc"),
						mdlwfdisc, am,
						chanmap,true);
				        }
        				/* We have to set up these tables because we are using wfprocess
        				to store 3c objects.*/
        				dbho.db.record=rec;
        				pwfid=dbho.get_int("pwfid");
        				rec=dbhoev.append();
        				dbhoev.put("pwfid",pwfid);
        				dbhoev.put("evid",evid);
        				rec=dbhosc.append();
        				dbhosc.put("sta",laststa);
        				/* Ignore the fact this has endian implications.  At present datatype
        				is used to handle this. */
        				dbhosc.put("chan","3C");
        				dbhosc.put("pwfid",pwfid);
        				/* This is necessary only if time is foobarred */
        				if(start_time_is_origin)
        				{
						/* These are frozen so will put them in */
						d3c.put("arrival.time",atime);
						d3c.put("arrival.iphase","P");
						d3c.put("assoc.phase","P");
						d3c.put("assoc.vmodel","iasp91");
						d3c.put("assoc.timedef","d");
						d3c.put("time",d3c.t0);
						d3c.put("endtime",d3c.endtime());
						/* The following should not be necessary, but ArrivalUpdater at this
						point does not handle aliases correctly */
						string sval;
						int ival;
						double dval;
						sval=d3c.get_string("sta");
						d3c.put("assoc.sta",sval);
						d3c.put("arrival.sta",sval);
						sval=d3c.get_string("chan");
						d3c.put("arrival.chan",sval);
						ival=d3c.get_int("orid");
						d3c.put("assoc.orid",ival);
						AUhandle.update(d3c);
                                        }
				}
				/* Cleanup and prep for new gather */
				tracecount=1;  // 1 because we push d_read below 
				ic=1; //ic is one for same reason as tracecount
				d.clear();
				d.push_back(d_read);
				laststa=sta;
				/* This must come from d_read because it is now current */
				lastevid=d_read.get_int("evid");
			}

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
