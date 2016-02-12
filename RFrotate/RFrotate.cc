#include <list>
#include <fstream>
#include "perf.h"
#include "coords.h"
#include "db.h"
#include "seispp.h"
#include "SphericalCoordinate.h"
#include "dbpp.h"
#include "Hypocenter.h"
#include "ThreeComponentSeismogram.h"
#include "ArrivalUpdater.h"
/*
Modification history: 
	Dec 15, 2014 by Xiaotao Yang
			1) printed out working progress for reference in verbose mode;
			2) added rotation method tag to wfprocess.algorithm appending RFrotate.
			3) commented out debugging lines for is_bundle test.
			4) debugged rotation method errors, should be lqt (and others) instead of -lqt (and -others).
*/

DatascopeHandle BuildCatalogView(DatascopeHandle& dbh)
{
	DatascopeHandle result(dbh);
	dbh.lookup("event");
        if(SEISPP_verbose)
            cerr << "Number of rows in event table="<<dbh.number_tuples()<<endl;
	list<string> sortkeys;
	sortkeys.push_back("evid");
	dbh.sort(sortkeys);
	dbh.natural_join("origin");
	string ss_to_prefor("orid==prefor");
	dbh.subset(ss_to_prefor);
        if(SEISPP_verbose)
            cerr << "Number of rows after join origin with (orid==prefor)="
                << dbh.number_tuples()<<endl;
	dbh.natural_join("assoc");
        if(SEISPP_verbose)
            cerr << "Number of rows after joint with assoc = "
                << dbh.number_tuples()<<endl;
	dbh.natural_join("arrival");
        if(SEISPP_verbose)
            cerr << "Final catalog view size (joining arrival)="
                << dbh.number_tuples()<<endl;
	return(dbh);
}
DatascopeHandle BuildWaveformView(DatascopeHandle& dbh)
{
    DatascopeHandle result(dbh);
    dbh.lookup("wfdisc");
    if(SEISPP_verbose)
            cerr << "BuildWaveformView:  "
                << "Number of rows in wfdisc="<<dbh.number_tuples()<<endl;
    dbh.natural_join("site");
    if(SEISPP_verbose)
            cerr << "BuildWaveformView:  "
                << "Number of rows after natural join of site table="
                <<dbh.number_tuples()<<endl;
    dbh.natural_join("sitechan");
    if(SEISPP_verbose)
            cerr << "BuildWaveformView:  "
                <<"Number of rows after joining sitechan="
                <<dbh.number_tuples()<<endl;
    return(dbh);
}
/* This builds working view by a join of the catalog view (dbcv)
   and the waveform view (dbwf) constructed using the functions
   above.  The view is then sorted and grouped for processing to
   build three-component seismograms. */
DatascopeHandle BuildWorkingView(DatascopeHandle& dbcv,
        DatascopeHandle& dbwv)
{
    try {
        list<string> jk1,jk2,sortkeys,groupkeys;
        DatascopeHandle result(dbwv);
        jk1.push_back(string("wfdisc.time::wfdisc.endtime"));
        jk1.push_back(string("sta"));
        jk2.push_back(string("arrival.time"));
        jk2.push_back(string("sta"));
        result.join(dbcv,jk1,jk2);
        if(SEISPP_verbose)
            cerr << "BuildWorkingView:  "
                << "Working view total final number of rows ="
                << result.number_tuples()<<endl;
        /*This is needed for 3c bundles */
        sortkeys.push_back(string("origin.time"));
        sortkeys.push_back(string("sta"));
        sortkeys.push_back(string("chan"));
        result.sort(sortkeys);
        cout << "Size of final view after sort="<<result.number_tuples()<<endl;
        
        //debug
        groupkeys.push_back(string("evid"));
        groupkeys.push_back(string("sta"));
        result.group(groupkeys);
        cout << "Size of final view after group="<<result.number_tuples()<<endl;
/*        if (result.is_bundle)
        	cout<<"DBhandle for working view is a bundle pointer."<<endl; */
        return(result);
    }
    catch(...){throw;};
}
Hypocenter extract_hypocenter(ThreeComponentSeismogram& d)
{
    double slat,slon,depth,otime;
    try {
        slat=d.get_double("origin.lat");
        slon=d.get_double("origin.lon");
        depth=d.get_double("origin.depth");
        otime=d.get_double("origin.time");
        /* Freeze method and model as tttaup and iasp91.  Not elegant
           but minor importance for now */
        return(Hypocenter(rad(slat),rad(slon),depth,otime,
                    string("tttaup"),string("iasp91")));
    }catch(SeisppError& serr)
    {
        throw serr;
    }
}
//test if the rms of the data passes the threshold set by the user.
bool rms_test(TimeSeries& d, double rms_max)
{
	double rms=dnrm2(d.s.size(),&(d.s[0]),1);
	if(rms>rms_max) 
	{
		d.live=false;
		return(false);
	}
	else
		return(true);
}
bool rms_test(ThreeComponentSeismogram& d, double rms_max)
{
	bool pass(true);
	for(int i=0; i<3; ++i)
	{
		TimeSeries * data=ExtractComponent(d,i);
		if(!rms_test(*data,rms_max))
		{
			d.live=false;
			pass=false;
		}
		delete data;
	}
	return(pass);
}     
void usage()
{
    cout<<"version <2.1> 8/19/2015"<<endl;
    cerr << "RFrotate nez2rtz|rtz2nez dbin dbout [-d outdir][-lqt|-fst -pf pffile  -v|V]"<<endl;
    cerr<<"    < nez2rtz > mode: -lqt or -fst methods could be used (default is -rtz)."<<endl
    	<<"    < rtz2nez > mode: for now, only data rotated using rtz method could be "<<endl
    	<<"                      rotated back correctly."<<endl;
    cerr<<"-d outdir: use alternate outdir, default is RFDataRotated."<<endl;
    cout<<endl<<"Xiaotao Yang & Gary Pavlis, Indiana University Bloomington"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    string pfname("RFrotate");
    /* parse the argument list*/
    if(argc<4) usage();
    string rotate_mode(argv[1]);
    if(rotate_mode!="nez2rtz" && rotate_mode!="rtz2nez")
    {
    	cerr<<"ERROR: Unknown rotate mode."<<endl;
    	usage();
    }
    string dbin_name(argv[2]);
    string dbout_name(argv[3]);
    string outdir("RFDataRotated");
    cout << "RFrotate:   Reading from database [ "
        << dbin_name<<" ]."<<endl
        << "	Writing results to database [ " <<dbout_name<<" ]."<< endl;
    int i,k;
    string rotation_method("rtz");
    
    for(i=4;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-lqt")
            rotation_method="lqt"; // replaced rotation_method=sarg. XT Yang
        else if(sarg=="-fst")
            rotation_method="fst"; // replaced rotation_method=sarg. XT Yang 
        else if(sarg=="-v" || sarg=="-V")
            SEISPP_verbose=true;
        else if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            sarg=string(argv[i]);
            pfname=sarg;
        }
        else if(sarg=="-d")
        {
            ++i;
            if(i>=argc) usage();
            sarg=string(argv[i]);
            outdir=sarg;
        }
        else
            usage();
    }
    /* Open a log file in append mode */
	ofstream logfile;
	logfile.open("RFrotate.log",ios::app);
	
    if(rotate_mode=="rtz2nez") rotation_method="rtz";
    if(SEISPP_verbose)
    {
    	cout<<"RFrotate mode: "<<rotate_mode<<endl;
		cout<<"rotation method: "<<rotation_method<<endl; // print out rotation method. XT Yang
	}
    /*Standard way to open Antelope pf file.*/
    Pf *pf;
    if(pfread(const_cast<char *>(pfname.c_str()),&pf))
    {
        cerr << "pfread failed on pf file named "<<pfname<<endl;
        exit(-1);
    }
    try {
        //PfStyleMetadata md=pfread(pfname);
        Metadata control(pf);
        /* We need this list of attributes to load from working
           view we pull from the parameter file.  */
        bool use_wfdisc_in=control.get_bool("use_wfdisc_in");
        MetadataList inmdl;
        if(use_wfdisc_in) 
        {
        	if(rotate_mode=="nez2rtz")
        		inmdl=pfget_mdlist(pf,string("input_mdl_wfdisc_nez2rtz"));
        	else
        		inmdl=pfget_mdlist(pf,string("input_mdl_wfdisc_rtz2nez"));
        }
        else
        {
        	if(rotate_mode=="nez2rtz")
        		inmdl=pfget_mdlist(pf,string("input_mdl_wfprocess_nez2rtz"));
        	else
        		inmdl=pfget_mdlist(pf,string("input_mdl_wfprocess_rtz2nez"));
        }
        MetadataList outmdl_wfd=pfget_mdlist(pf,string("output_wfdisc_mdl"));
        MetadataList outmdl_wfp=pfget_mdlist(pf,string("output_wfprocess_mdl"));
        MetadataList mdlassoc;
		MetadataList mdlarrival;
		
		//parameters used in rtz2nez.
		bool start_time_is_origin;
		string rchan("-");
		string tchan("-");
		string zchan("-");
		int nchan_max(2);
		bool use_vertical_chan(false);
		/* QC problem with EARS data found a problem with high 
		   amplitude traces that dominated the stack. This parameter
		   will autokill (with a message) any trace with an rms amplitude
		   exceeding this size on either radial or transverse */
		double rms_kill_threshold;
		if(rotate_mode=="rtz2nez")
		{
			mdlassoc=pfget_mdlist(pf,"save_assoc_mdl");
			mdlarrival=pfget_mdlist(pf,"save_arrival_mdl");
			start_time_is_origin=control.get_bool("start_time_is_origin");
			rchan=control.get_string("radial_channel_code");
			tchan=control.get_string("transverse_channel_code");
			use_vertical_chan=control.get_bool("use_vertical_chan");
			if(use_vertical_chan)
			{	zchan=control.get_string("vertical_channel_code"); nchan_max=3;}
			//save_wfdisc=control.get_bool("save_wfdisc");
			rms_kill_threshold=control.get_double("rms_kill_threshold");
		}
        vector<string> output_channels;
        output_channels.push_back(control.get_string("out_chan_code_x1"));
        output_channels.push_back(control.get_string("out_chan_code_x2"));
        output_channels.push_back(control.get_string("out_chan_code_x3"));
        AttributeMap am("css3.0");
        bool save_wfdisc=control.get_bool("save_wfdisc_table");
        bool save_wfprocess=control.get_bool("save_wfprocess_table");
        if(!save_wfdisc && !save_wfprocess)
        {
        	cerr<<"ERROR: you have to save at least one of the wfdisc and wfprocess tables."<<endl;
        	exit(-1);
        }
        //vp0, vs0 - P and S velocities of the half space = surface
        double vp0,vs0;
        vp0=control.get_double("vp0");
        vs0=control.get_double("vs0");
        //string outdir=control.get_string("OutPutDirectory");
        string outdfilebase=control.get_string("OutPutFileNameBase");

        /* First prepare all the database handles */
        DatascopeHandle dbin(dbin_name,true);
        DatascopeHandle dbout(dbout_name,false);
        
        //debug
        //DatascopeHandle dbotmp("dbouttmp",false);
        if(rotate_mode=="nez2rtz")
        {
			dbin=BuildCatalogView(dbin);
			DatascopeHandle dbwf(dbin);
			dbwf=BuildWaveformView(dbwf);
			dbin=BuildWorkingView(dbin,dbwf);
			dbout.lookup("wfdisc");
			DatascopeHandle dbo3c(dbout);
			dbo3c.lookup("wfprocess");
			DatascopeHandle dboevl(dbout);
			dboevl.lookup("evlink");
			DatascopeHandle dboscl(dbout);
			dboscl.lookup("sclink");
		
			// print out progress.  XT Yang
			cout<<"------------------------------"<<endl
				<<"Starting rotation process ..."<<endl;
			/* We loop over the group pointer to get one 3c seismogram 
			   at a time */
			dbin.rewind();

			int nrows=dbin.number_tuples();
			string sta;
			for(i=0;i<nrows;++i,++dbin)
			//for(i=0;i<1;++i,++dbin) // debug
			{
				try
				{
					/*	if (dbin.is_bundle)
						cout<<"DBhandle for dbin is a bundle pointer."<<endl; */

					// print out progress.  XT Yang
					cout<<"Working on row [ "<<i+1<<" ] of [ "<<nrows<<" ]."<<endl;
					ThreeComponentSeismogram d(dbin,inmdl,am); 
					/* This may not be essential, but better safe than sorry */
					d.rotate_to_standard();
			
					/* We define this here so we can use it in messages later
					   if needed and for sclink below*/
					sta=d.get_string("sta");
					Hypocenter h(extract_hypocenter(d));
					double rlat,rlon,relev;
					rlat=d.get_double("site.lat");
					rlon=d.get_double("site.lon");
					relev=d.get_double("site.elev");
					rlat=rad(rlat);  rlon=rad(rlon);
					/* For now assume P wave phase for slowness vector when
					   the fst is used */
					SlownessVector u=h.pslow(rlat,rlon,relev);
					//DEBUGGING
					//cout<<"slowness: ux, uy ="<<u.ux<<", "<<u.uy<<endl;
			
					double az=u.azimuth();
					SphericalCoordinate sc_for_rotation;
					/*
					The data are rotated such that x1 becomes 
		 				the transverse component, x2 becomes radial, 
		 				and x3 becomes longitudinal.  
		 			Commented by Xiaotao Yang.
					*/
			
					if(rotation_method=="fst")
					{
						d.free_surface_transformation(u,vp0,vs0);
					}
					else if(rotation_method=="lqt")
					{
						sc_for_rotation=PMHalfspaceModel(vp0,vs0,u.ux,u.uy);
						//DEBUGGING
					/*cout<<"station = "<<sta<<endl;
					//cout<<"wfdisc time = "<<d.get_
					cout<<"phi ="<<sc_for_rotation.phi<<endl;
					cout<<"theta ="<<sc_for_rotation.theta<<endl;
					cout<<"M_PI_2 - phi = "<<M_PI_2-sc_for_rotation.phi<<endl;*/
						d.rotate(sc_for_rotation);
					}
					else
					{
						sc_for_rotation.radius=1.0;
						sc_for_rotation.theta=0.0;
						sc_for_rotation.phi=M_PI_2-az;
						d.rotate(sc_for_rotation);
					}
					/* Now we unbundle the data and write it to output db*/
					d.put("dir",outdir);
					string sdfile;
					sdfile=outdfilebase+"_"+sta+".3C";
					d.put("dfile",sdfile);
					if(save_wfprocess)					
					{	// add rotation method to wfprocess.algorithm. XT Yang
						string wp_algorithm="RFrotate("+rotation_method+")";
						d.put("wfprocess.algorithm",wp_algorithm);
						//    d.put("wfprocess.algorithm",string("RFrotate"));
						d.put("timetype",string("a"));
						/* Not certain this is necessary, but minimal cost */
						//d.put("dir",dir3c);
						//d.put("dfile",dfile3c);
						long pwfid;
						pwfid=dbsave_oriented(d,dbo3c.db,string("wfprocess"),outmdl_wfp,am);
						/* We trap these errors and abort because if this section fails
						   there is a setup problem. */
						int retcode;
						long evid=d.get_long("evid");
						retcode=dbaddv(dboevl.db,0,"evid",evid,"pwfid",pwfid,NULL);
						if(retcode==dbINVALID)
						{
							cerr << "Fatal:  dbaddv failed writing evlink table"<<endl
								<< "evid="<<evid<<" pwfid="<<pwfid<<endl
								<< "Check parameter file and required database tables"
								<<endl;
							exit(-1);
						}
						const string chanout("BUNDLE3C");
						retcode=dbaddv(dboscl.db,0,"sta",sta.c_str(),
								"chan",chanout.c_str(),
								"pwfid",pwfid,NULL);
						if(retcode==dbINVALID)
						{
							cerr << "Fatal:  dbaddv failed writing sclink table"<<endl
								<< "sta="<<sta
								<< " evid="<<evid<<" pwfid="<<pwfid<<endl
								<< "Check parameter file and required database tables"
								<<endl;
							exit(-1);
						} 
					} 
					if(save_wfdisc)
					{
						TimeSeries *component;
						vector<string>::iterator iptr;
						for(k=0,iptr=output_channels.begin();k<3;k++,iptr++) 
						{	
							//debug
							//cout<<"comp="<<k<<", chan="<<*iptr<<endl;
							component=ExtractComponent(d,k);
							component->put("chan",*iptr);
							dbsave(*component,dbout.db,string("wfdisc"),outmdl_wfd,am);
							delete component;
						}
					}         
					//print out progress. XT Yang
					if(SEISPP_verbose)
						cout<<"Done processing row [ "<<i+1<<" ] of [ "<<nrows<<" ]."<<endl;
				}catch(SeisppError& serr)
				{
				  cerr << "Something threw a SeisppError exception in main "
					  <<"processing loop for station="<<sta
					  <<endl<<"Error message follows:"<<endl;
				  serr.log_error();
				}
			} 
            
        }//end of nez2rtz.
        else  //start of rtz2nez rotation process.
		{
			try 
			{
				DatascopeHandle dbhi(dbin_name,true);
				if(use_wfdisc_in)
				{
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
					sortkeys.push_back("sta");
					sortkeys.push_back("origin.time");
					dbhi.sort(sortkeys);
				}
				else
				{
					dbhi.lookup("evlink");
					dbhi.natural_join("origin");
					dbhi.natural_join("sclink");
					dbhi.natural_join("wfprocess");
					list<string> jk1, jk2;
					jk1.push_back(string("sta"));
					jk2.push_back(string("sta"));
					dbhi.leftjoin("site",jk2,jk1);
					list<string> sortkeys;
					sortkeys.push_back("sta");
					sortkeys.push_back("origin.time");
					dbhi.sort(sortkeys);
				}

				cout << "Working input view has "<<dbhi.number_tuples()
					<<" rows"<<endl;
				//DEBUG
				//exit(-1);
				if(dbhi.number_tuples()<=0)
				{
					if(use_wfdisc_in)
						cerr << "No data to process.  Exiting.  "<<endl
							<< "Working view "
							<<"(wfdisc<-event->origin->assoc->arrival->site)"
							<< " is empty"<<endl;
					else
						cerr << "No data to process.  Exiting.  "<<endl
							<< "Working view "
							<<"(evlink->origin->sclink->wfprocess->site)"
							<< " is empty"<<endl;
					exit(-1);
				}
				DatascopeHandle dbho(dbout_name,false);
				dbho.lookup("wfprocess");
				/* Need to put stuff in evlink and sclink for pwmig */
				DatascopeHandle dbhosc(dbho);
				dbhosc.lookup("sclink");
				DatascopeHandle dbhoev(dbho);
				dbhoev.lookup("evlink");
				DatascopeHandle dbhwfdisc(dbho);
				dbhwfdisc.lookup("wfdisc");
				/* We use this as a convenient way to create a new assoc and arrival */
				/*
				//not deal with this issue for now. Xiaotao Yang
				ArrivalUpdater AUhandle(dynamic_cast<DatabaseHandle&>(dbho),
						mdlassoc,mdlarrival,string("css3.0"));
				*/ 
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
				//string outdir("wf3c");
				//string dfilebase("RFdata");
				dbhi.rewind();
				laststa=dbhi.get_string("sta");
				lastevid=dbhi.get_int("evid");
				/* This holds radial and transverse - hence 2*/
				vector<TimeSeries> d;
				int ic(0);  // always 0 or 1 */
				int tracecount(0);
				int ntuples=dbhi.number_tuples();
				bool datatype3c(false);
				string datatype=dbhi.get_string("datatype");
				if(datatype=="c3") datatype3c=true;
				
				for(i=0, dbhi.rewind();i<ntuples;++i,++dbhi)
				{
					bool datalive(true);
					bool ready_to_save(false);
					ThreeComponentSeismogram d3c;
					cerr<<"----------------- [ "<<i+1<<" / "<<ntuples<<" ]"<<endl;
					if(!datatype3c)  //for 1c data (time series scalar data).
					{
						/* This holds current time series data read.  Accumulate until sta or event change */
						//cerr<<"flag0"<<endl;
						TimeSeries d_read(dynamic_cast<DatabaseHandle&>(dbhi),inmdl,am);
						//cerr<<"flag 1"<<endl;
						/* ic must index member of d just pushed. Loop assumes ic=0 on first pass
						and is incremented thereafter.  Reset to 0 when new event-sta encounterd */
						sta=d_read.get_string("sta");
						evid=d_read.get_int("evid");
						//cerr<<"flag 2"<<endl;
						/* We compute rms here assuming we have no 
						   gaps. There is a fancier way to do this in
						   the seispp library to deal with gaps but that
						   baggage is not needed here */
						if(SEISPP_verbose) 
							cout << "Indices i,sta,chan,evid:  "
								<<i<<", "<<sta
								<<", "<<d_read.get_string("chan")
								<<", "<<evid
								<<"."<<endl;
						if(!rms_test(d_read,rms_kill_threshold))
						{
							cout << "High rms for "<<sta<<" on channel "
								<<d_read.get_string("chan")<<endl
								<< "Marking for deletion"<<endl;
							d_read.live=false;
						}
						if(sta==laststa && evid==lastevid && (i != (ntuples-1)))
						{
							if(tracecount>nchan_max)
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
							slat=d[ic].get_double("origin.lat"); // ic always 1 in this section
							slon=d[ic].get_double("origin.lon");
							evid=d[ic].get_int("evid");
							cout << "Processing event "<<evid
								<< " located at (lat,lon)=("<<slat
								<<", "<<slon<<") "
								<<" for station="<<laststa<<endl;
	
							slat=rad(slat);
							slon=rad(slon);
							sdepth=d[ic].get_double("origin.depth");
							stime=d[ic].get_double("origin.time");
							Hypocenter h(slat,slon,sdepth,stime,method,model);
							stalat=d[ic].get_double("site.lat");
							stalon=d[ic].get_double("site.lon");
							staelev=d[ic].get_double("site.elev");
							stalat=rad(stalat);
							stalon=rad(stalon);
							
							if(use_wfdisc_in)
							{
								atime=d[ic].get_double("arrival.time");
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
							azimuth.  In seispp the transformation is the forward transformation so
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
							ThreeComponentSeismogram d3c0(dynamic_cast<Metadata&>(d[ic]),false);
							d3c0.t0=d[ic].t0;
							d3c0.ns=d[ic].ns;
							d3c0.dt=d[ic].dt;
							d3c0.components_are_orthogonal = true;
							/* initialize */
							for(k=0;k<d[ic].ns;++k)
								for(j=0;j<3;++j) d3c0.u(j,k)=0.0;
							/* stuff r into 1 (x2) and t into 0 (x1).  Leave x3 0 */
							/*
							Commented by Xiaotao Yang
							added option to handle three components, i.e., use vertical.
							*/
							int iic,jc;
						
							for(iic=0,datalive=true;iic<tracecount;++iic)
							{
								if(d[iic].live)
								{
									string chan;
									chan=d[iic].get_string("chan");
									if(chan.find(rchan)!=string::npos)
										jc=1;
									else if(chan.find(tchan)!=string::npos)
										jc=0;
									else if(chan.find(zchan)!=string::npos)
										jc=2;
									else
									{
										cerr << "Fatal: do not know how to handle channel "
											<< chan <<endl;
										exit(-1);
									}
									for(k=0;k<d[iic].ns;++k) d3c0.u(jc,k)=d[iic].s[k];
								}
								else
								{
									datalive=false;
									break;
								}
							}
							d3c=d3c0;
							ready_to_save=true;
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
					else  //for 3c datatype.
					{
						ThreeComponentSeismogram d3c0(dbhi,inmdl,am);
						//DEBUG
						/*
						cout << dynamic_cast<Metadata&>(d3c0);
						dmatrix tmp=d3c0.u.tr();
						cout << tmp<<endl;
						*/
						sta=d3c0.get_string("sta");
						evid=d3c0.get_int("evid");
						cerr<<"Working on station = "<<sta
							<<", for evid = "<<evid<<endl;
						
						if(!rms_test(d3c0,rms_kill_threshold))
						{
							cerr<<"* Data has high rms. Marked for deletion."<<endl;
							d3c0.live=false;
							datalive=false;
						}
						
						slat=d3c0.get_double("origin.lat"); // ic always 1 in this section
						slon=d3c0.get_double("origin.lon");
						evid=d3c0.get_int("evid");
						cout<< " located at (lat,lon)=("<<slat
							<<", "<<slon<<") "
							<<" for station="<<sta<<endl;

						slat=rad(slat);
						slon=rad(slon);
						sdepth=d3c0.get_double("origin.depth");
						stime=d3c0.get_double("origin.time");
						Hypocenter h(slat,slon,sdepth,stime,method,model);
						stalat=d3c0.get_double("site.lat");
						stalon=d3c0.get_double("site.lon");
						staelev=d3c0.get_double("site.elev");
						stalat=rad(stalat);
						stalon=rad(stalon);
						
						seaz=h.seaz(stalat,stalon);
						logfile<<"evid = "<<evid<<", sta = "<<sta<<", seaz = "<<180*seaz/M_PI<<endl;
						/* radial is defined as propagation direction, but seaz is back azimuth */
						az=seaz-M_PI;
						if(az<(-2.0*M_PI)) az+=(2.0*M_PI);
						if(SEISPP_verbose) cout <<" station="<<sta
							<<" computed radial azimuth="<<deg(az)
							<<" Distance(deg)="<<deg(h.distance(stalat,stalon))<<endl;
						/* Now we have to set the transformation matrix from az.  Convention here
						is we put R into x2, so the rotation angle turns out to be -az when you 
						work through the difference between rotation from the x1 axis and an
						azimuth.  In seispp the transformation is the forward transformation so
						the signs for a and b are as list */
						a=cos(-az);
						b=sin(-az);
						//DEBUG
						/*
						cerr<<"Before put: "<<endl;
						cerr<<"U11 = "<<d3c0.get_double("U11")<<endl
							<<"U21 = "<<d3c0.get_double("U21")<<endl
						<<"U31 = "<<d3c0.get_double("U31")<<endl
						<<"U12 = "<<d3c0.get_double("U12")<<endl
						<<"U22 = "<<d3c0.get_double("U22")<<endl
						<<"U32 = "<<d3c0.get_double("U32")<<endl
						<<"U13 = "<<d3c0.get_double("U13")<<endl
						<<"U23 = "<<d3c0.get_double("U23")<<endl
						<<"U33 = "<<d3c0.get_double("U33")<<endl;
						*/
						/* by loading these in the metadata for d they will be cloned
						and used to construct the transformation matrix for d3c. */
						////
						d3c0.put("U11",a);
						d3c0.put("U21",-b);
						d3c0.put("U31",0.0);
						d3c0.put("U12",b);
						d3c0.put("U22",a);
						d3c0.put("U32",0.0);
						d3c0.put("U13",0.0);
						d3c0.put("U23",0.0);
						d3c0.put("U33",1.0);
						
						//DEBUG
						/*
						cerr<<"What to put: "<<endl;
						cerr<<"U11 = "<<a<<endl
							<<"U21 = "<<-b<<endl
						<<"U31 = "<<0.0<<endl
						<<"U12 = "<<b<<endl
						<<"U22 = "<<a<<endl
						<<"U32 = "<<0.0<<endl
						<<"U13 = "<<0.0<<endl
						<<"U23 = "<<0.0<<endl
						<<"U33 = "<<1.0<<endl;
						*/
						
						ThreeComponentSeismogram d3c_tmp(dynamic_cast<Metadata&>(d3c0),false);
						d3c_tmp.u=d3c0.u;
						d3c_tmp.components_are_orthogonal = true;
						d3c=d3c_tmp;
						//debug
						//cerr<<"tmatrix 00 after put: "<<d3c.tmatrix[0][0]<<endl;
						ready_to_save=true;
						datalive=true;
						//DEBUG
						/*
						cerr<<"After (retrieved) put: "<<endl;
						cerr<<"U11 = "<<d3c.get_double("U11")<<endl
							<<"U21 = "<<d3c.get_double("U21")<<endl
						<<"U31 = "<<d3c.get_double("U31")<<endl
						<<"U12 = "<<d3c.get_double("U12")<<endl
						<<"U22 = "<<d3c.get_double("U22")<<endl
						<<"U32 = "<<d3c.get_double("U32")<<endl
						<<"U13 = "<<d3c.get_double("U13")<<endl
						<<"U23 = "<<d3c.get_double("U23")<<endl
						<<"U33 = "<<d3c.get_double("U33")<<endl;
						*/
						//DEBUG
						//exit(-1);
					}
					
					//save to db
					if(datalive && ready_to_save) 
					{
						d3c.live=true;
						d3c.rotate_to_standard();
						/* We need to post these things or we're screwed */
						d3c.put("dir",outdir);
						char buf[128];
						string ss=string(outdfilebase)+"_"+string(sta)+".w3c";
						d3c.put("dfile",string(ss));
						
						d3c.put("wfprocess.algorithm",(char *)"RFR_rtz2nez");
						d3c.put("timetype",(char *)"a");
						//cerr<<"test save"<<endl;
						//cout << "3c data for sta "<<d3c.get_string("sta")<<" t0="<<strtime(d3c.t0)<<endl;
						if(save_wfprocess)
						{
							int rec=dbsave(d3c,dbho.db,string("wfprocess"),outmdl_wfp,am);
							if(rec<0)
							{
								cerr << "dbsave failed working on input record="<<i<<endl;
								exit(-1);
							}
							/* We have to set up these tables because we are using wfprocess
							to store 3c objects.*/
							dbho.db.record=rec;
							pwfid=dbho.get_int("pwfid");
							rec=dbhoev.append();
							dbhoev.put("pwfid",pwfid);
							dbhoev.put("evid",evid);
							rec=dbhosc.append();
							dbhosc.put("sta",sta);
							/* Ignore the fact this has endian implications.  At present datatype
							is used to handle this. */
							dbhosc.put("chan","3C");
							dbhosc.put("pwfid",pwfid);
						}
						
						if(save_wfdisc)
						{
							ss=string(outdfilebase)+"_"+string(sta)+".w";
							d3c.put("dfile",string(ss));
							dbsave(d3c,dbhwfdisc.db,
							string("wfdisc"),
							outmdl_wfd, am,
							output_channels,true);
						}
						/* This is necessary only if time is foobarred */
						/*
						if(start_time_is_origin && use_wfdisc_in)
						{
							// These are frozen so will put them in
							d3c.put("arrival.time",atime);
							d3c.put("arrival.iphase",(char *)"P");
							d3c.put("assoc.phase",(char *)"P");
							d3c.put("assoc.vmodel",(char *)"iasp91");
							d3c.put("assoc.timedef",(char *)"d");
							d3c.put("time",d3c.t0);
							d3c.put("endtime",d3c.endtime());
							// The following should not be necessary, but ArrivalUpdater at this
							//point does not handle aliases correctly
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
						*/
					}
				}
			}catch (SeisppError& serr)
			{
				serr.log_error();
			}
		}
		
    }catch(SeisppError& serr)
    {
        cerr << "Fatal error before processing started"<<endl;
        serr.log_error();
    }
    logfile.close();
}
