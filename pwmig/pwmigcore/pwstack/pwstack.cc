#include <string>
#include <sstream>
#include "stock.h"
//#include "glputil.h"
#include "gclgrid.h"
#include "Hypocenter.h"
#include "seispp.h"
#include "PwmigFileHandle.h"
#include "PfStyleMetadata.h"
#include "SlownessVectorMatrix.h"
#include "pwstack.h"
#include "pwstack_reader.h"

using namespace std;
using namespace SEISPP;

bool Verbose;

void usage()
{
    cerr << "Usage:  pwstack datafile [-np np -rank rank -v -pf pfname]"
        <<endl <<"  If -np is set every np events processed starting at rank"
        <<endl <<"  -v - run verbose"
        <<endl <<"  -pf use alternate parameter file pfname"
        <<endl;
    exit(-1);
}


void load_file_globals(PwmigFileHandle& fh,int evid, double olat, double olon, 
	double odepth, double otime, string gridname)
{
	fh.filehdr.evid=evid;
	fh.filehdr.slat=olat;
	fh.filehdr.slon=olon;
	fh.filehdr.sdepth=odepth;
	fh.filehdr.stime=otime;
	strncpy(fh.filehdr.gridname,gridname.c_str(),16);
}
/* Removes members of d that have no data inside time aligned window range */
auto_ptr<ThreeComponentEnsemble> clean_gather(
	auto_ptr<ThreeComponentEnsemble> d,
		double tstart, double tend)
{
	vector<ThreeComponentSeismogram>::iterator dptr;
	for(dptr=d->member.begin();dptr!=d->member.end();++dptr)
	{
		//DEGUG
		/*
		cout << "clean_gather sta="<<dptr->get_string("sta")
			<< " t0="<< (dptr->t0) <<endl;
		*/
		double t0d=dptr->t0;
		double ted=dptr->endtime();
		/* Test for data outside stack window and delete
 		them.  dptr assignment construct is used because
		stl erase method returns iterator of next element
		after the delete or the end.  Proper way to 
		do this with an stl container. */
		if((t0d>tend)||(ted<tstart))
		{
			if(SEISPP_verbose) cerr << "clean_gather:  "
					<< "deleting seismogram for sta="
					<< dptr->get_string("sta")
					<<endl;
			dptr=d->member.erase(dptr);
		}
	}
	return(d);
}

	



bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
    int i,j;
    char *pfin=NULL;
    char *ensemble_tag;

    // Tbl tag used in pf to define depth-dependent apeture.
    // frozen here as this constant but passed as a ariable to the
    // appropriate constructor below
    string aperture_tag="depth_dependent_aperture";
    string mdlist_tag="copy_metadata_list";

    ios::sync_with_stdio();
    /* Initialize the error log and write a version notice */

    /* usual cracking of command line */
    if(argc < 2) usage();
    string infile(argv[1]);
    int rank(0);
    int np(1);

    for(i=2;i<argc;++i)
    {
	string sarg(argv[i]);
	if(sarg=="-v")
        {
            Verbose=true;
            SEISPP::SEISPP_verbose=true;
        }
	else if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pfin = argv[i];
        }
        else if(sarg=="-rank")
        {
            ++i;
            if(i>=argc) usage();
            rank=atoi(argv[i]);
	    if(rank<0) 
	    {
		cerr << "Illegal value for -r argument ="
			<< rank<<endl
			<< "Must be nonnegative integer"<<endl;
		usage();
	    }
        }
	else if(sarg=="-np")
	{
	    ++i;
	    if(i>=argc) usage();
	    np=atoi(argv[i]);
	    if(np<=0) 
	    {
		cerr << "Illegal value for -np argument ="
			<< np<<endl
			<< "Must be positive integer"<<endl;
		usage();
            }
	}
        else
            usage();
    }
	/* sanity checks on manual parallel variables */
	if( (np>0) || (rank>=0) )
	{
		if( (np<0) || (rank<0) )
		{
			cerr << "Illegal argument combination.  "
				<< "-rank and -np must both be set if either is set"
				<< endl;
			usage();
		}
		else if(rank>(np-1))
		{
			cerr << "Illegal argument combination.  "
				<< "Value for -rank must be >= -np"<<endl;
			usage();
		}
	}
	else
	{
		/* single processor mode set if we end up here.*/
		np=1;
		rank=0;
	}
    /* this sets defaults */
    if(pfin == NULL) pfin = strdup("pwstack.pf");

    try
    {
	PfStyleMetadata control=SEISPP::pfread(string(pfin));
        // This builds the grid of plane wave components
        RectangularSlownessGrid ugrid(control,"Slowness_Grid_Definition");
        // control parameters on stack process
        double ts,te;
	ts=control.get_double("stack_time_start");
	te=control.get_double("stack_time_end");
	int stack_count_cutoff=control.get_int("stack_count_cutoff");
	double aperture_taper_length=control.get_double("aperture_taper_length");
	double centroid_cutoff=control.get_double("centroid_cutoff");
        // the data are windowed around arrivals to this interval
        // normally should be ts:te with sufficient allowance for moveout
        double tsfull, tefull;                    // normally longer than ts and te to allow
        tsfull = control.get_double("data_time_window_start");
        tefull = control.get_double("data_time_window_end");
        TimeWindow data_window(tsfull,tefull);
        string dir = control.get_string("waveform_directory");
        if(dir.size()<=0)
            dir=strdup("./pwstack");
        if(makedir(const_cast<char *>(dir.c_str())))
	{
	    cerr << "Cannot create directory ="<<dir<<endl;
            exit(-1);
	}
	string dfilebase=control.get_string("output_data_file_base_name");
	string fielddir=control.get_string("field_file_directory");
	string fieldnamebase=control.get_string("field_base_name");
        //
        // station_mdl defines database attributes copied to
        // metadata space of data object.  ensemble_mdl are
        // globals copied to the metadata area for the entire ensemble
        //
        MetadataList station_mdl=get_mdlist(control,"station_metadata");
        MetadataList ensemble_mdl=get_mdlist(control,"ensemble_metadata");
        bool use_fresnel_aperture=control.get_bool("use_fresnel_aperture");
        DepthDependentAperture aperture;
        if(use_fresnel_aperture)
        {
            cout << "Using Fresnel zone aperture widths for stacks"<<endl;
            double fvs,fvp,fcm,fperiod,fdt;
            int fnt;
            fvs=control.get_double("fresnel_vs");
            fvp=control.get_double("fresnel_vp");
            fcm=control.get_double("fresnel_cutoff_multiplier");
            fperiod=control.get_double("fresnel_period");
            fdt=control.get_double("fresnel_lag_time_sampling_interval");
            fnt=control.get_int("fresnel_number_lag_samples");
            cout << "Fresnel zone parameters:  Vs="<<fvs
                << ", Vp="<<fvp
                << ", period="<<fperiod
                << ", cutoff multiplier="<<fcm
                << "Aperture parameters will be specified by "<<fnt
                << " points sampled at intervals of " <<fdt
                << " for total depth of "<< fdt*static_cast<double>(fnt)
                <<endl;
            aperture=DepthDependentAperture(fvs,fvp,fperiod,fdt,fnt,
                    fcm,true);
        }
        else 
        {
            aperture=DepthDependentAperture(control,aperture_tag);
        }
	/* These will throw an exception if not defined.  
	We use a parameter to turn these off if desired.
	Not very elegant, but functional. */
        TopMute mute(control,string("Data_Top_Mute"));
	bool enable_data_mute=control.get_bool("enable_data_mute");
	if(enable_data_mute)
		mute.enabled=true;
	else
		mute.enabled=false;
        TopMute stackmute(control,string("Stack_Top_Mute"));
	bool enable_stack_mute=control.get_bool("enable_stack_mute");
	if(enable_stack_mute)
		stackmute.enabled=true;
	else
		stackmute.enabled=false;
        /* new parameters April 2007 for coherence attribute calculator*/
        double dtcoh,cohwinlen;
        dtcoh=control.get_double("coherence_sample_interval");
        cohwinlen=control.get_double("coherence_average_length");
        /* This access a special format input data file for reading.*/
        PwstackBinaryFileReader input_handle(infile);
	cout << "Processing begins on file " 
		<<  infile << endl
		<<"Number of events in file= "
                <<input_handle.number_events()
                <<endl;

        // We need to load the primary GCLgrid that defines the
        // location of pseudostation points.  It is assumed we
        // will build a site and sitechan table in the output database
        // in that define these virtual stations in a different
        // application.  i.e. dbho will eventually need a site table
        // produced separately to be useful.  This program should not
        // need this information though.
        //
        string stagridname=control.get_string("pseudostation_grid_name");
        /* This is a file based constructor with two implied names
        ending in .dat and .hdr */
        GCLgrid stagrid(stagridname);
        cout << "Using pseudostation grid of size "<<stagrid.n1<<" X " <<stagrid.n2<<endl;

	/* This field derived from stagrid is used to hold stack fold
	for each pseudostation grid point */
	GCLscalarfield fold(stagrid);
	fold.zero();  // best to guarantee initialization

	char ssbuf[256];
        int event_number;
        for(event_number=rank;event_number<np;event_number+=np)
        {
            int iret;
            long evid;
            double olat,olon,odepth,otime;
            // ensemble is read once for entire grid in this
            // version of the program.  This assumes newer
            // methods like that under development by Fan
            // will prove better than pseudostation method
            //
            auto_ptr<ThreeComponentEnsemble> ensemble;
            try{
              ensemble=auto_ptr<ThreeComponentEnsemble>(input_handle.read_gather(event_number));
            }
            catch(std::exception& stdexc)
            {
                cerr << "Error reading ensemble for evid="<<evid<<endl
                    << "Message from reader:"<<endl
                    << stdexc.what()<<endl
                    << "Skipping data for this event"<<endl;
                continue;
            }
	    evid=ensemble->get_long("evid");
            /* This was added 2015 to remove dependence on global travel
               time calculators.  Previously we computed slowness vectors
               in this program.  Now we import them through this mechanism */
            SlownessVectorMatrix svm=input_handle.gather_svm();
            /* Chaos will result if we don't require the svm to match
               the dimensions of the stagrid.  svm has no coordinates so
               the only sanity check is dimensions */
            if( (svm.rows() != stagrid.n1) || (svm.columns()!=stagrid.n2) )
            {
                cerr << "WARNING:  semifatal error for ensemble with evid="
                    <<evid<<endl
                    << "Size mismatch between pseudostation gridname= "
                    << stagridname<<endl
                    << stagridname << " dimensions are "<< stagrid.n1
                    <<"X"<<stagrid.n2<<endl
                    <<"SlownessVectorMatrix in input data file dimension is "
                    <<svm.rows()<<"X"<<svm.columns()<<endl
                    <<"This ensemble was not processed"<<endl;
                continue;
            }


            /* Used to need this - retain as reminder until debug finished
            auto_ptr<ThreeComponentEnsemble>
                ensemble=ArrivalTimeReference(*din,"arrival.time",
                data_window);
                */
            // Release this potentially large memory area
	    if(SEISPP_verbose) cout << "Ensemble for evid="<<evid
			<< " has "<<ensemble->member.size()<<" seismograms"
			<<endl;
            /* Throws out data without overlaping times*/
	    int original_ensemble_size=ensemble->member.size();
	    ensemble=clean_gather(ensemble,ts,te);
	    int clean_ensemble_size=ensemble->member.size();
	    if(clean_ensemble_size!=original_ensemble_size)
	    {
		cerr << "Warning:  potential data problem"<<endl
			<< "clean_gather procedure cleared "
			<<(clean_ensemble_size-original_ensemble_size)
			<<" seismograms with time range outside "
			<<"stack window"<<endl
			<<"Number of seismograms remaining="
			<<ensemble->member.size()<<endl;
	    }
            // this should probably be in a try block, but we need to
            // extract it here or we extract it many times later.
            evid=ensemble->get_int("evid");
            olat=ensemble->get_double("lat");
            olon=ensemble->get_double("lon");
            // Database has these in degrees, but we need them in radians here.
            olat=rad(olat);  olon=rad(olon);
            odepth=ensemble->get_double("depth");
            otime=ensemble->get_double("time");
	    /* A way to assure this is cleared.  May not be necessary */
	    ssbuf[0]='\0';
	    stringstream ss(ssbuf);
	    ss << dir << "/" << dfilebase << "_" << evid;
	    string dfilebase=ss.str();
	    string dfile=dfilebase+DataFileExtension;
	    string coh3cf=dfilebase+Coh3CExtension;
	    string cohf=dfilebase+CohExtension;
	    PwmigFileHandle dfh(dfile,false,ugrid);
	    load_file_globals(dfh,evid,olat,olon,odepth,otime,stagridname);
	    PwmigFileHandle coh3cfh(coh3cf,false,ugrid);
	    load_file_globals(coh3cfh,evid,olat,olon,odepth,otime,stagridname);
	    PwmigFileHandle cohfh(cohf,true,ugrid);
	    load_file_globals(cohfh,evid,olat,olon,odepth,otime,stagridname);
            double lat0,lon0,elev0;
	    cout << "Beginning processing for event id = "<<evid<<endl;
            for(i=0;i<stagrid.n1;++i)
                for(j=0;j<stagrid.n2;++j)
            {
                lat0=stagrid.lat(i,j);
                lon0=stagrid.lon(i,j);
                elev0=stagrid.depth(i,j);
                elev0=-elev0;                     // z is positive down, elev is positive up
                // ux0 and uy0 are the incident wavefield's
                // slowness vector.  Stacks are made relative to this
                // and the data are aligned using P wave arrival times
                try
                {
                    //
                    // this is a backdoor method to pass this information
                    // to the main processing program called just below.
                    // perhaps bad form, but it is logical as these
                    // are metadata related to these data by any measure.
                    ensemble->put("ix1",i);
                    ensemble->put("ix2",j);
                    ensemble->put("lat0",lat0);
                    ensemble->put("lon0",lon0);
                    ensemble->put("elev0",elev0);
                    ensemble->put("gridname",stagridname);
                    /*
                    Hypocenter hypo(olat,olon,odepth,otime,
                        string("tttaup"),string("iasp91"));
                    SlownessVector slow=hypo.pslow(lat0,lon0,elev0);
                        */
                    SlownessVector slow=svm(i,j);
                    ensemble->put("ux0",slow.ux);
                    ensemble->put("uy0",slow.uy);

                    iret=pwstack_ensemble(*ensemble,
                        ugrid,
                        mute,
                        stackmute,
			stack_count_cutoff,
                        ts,
                        te,
                        aperture,
			aperture_taper_length,
			centroid_cutoff,
                        dtcoh,
                        cohwinlen,
                        ensemble_mdl,
                        dfh,
			coh3cfh,
			cohfh);
		    if(iret<=0)
		    {
			fold.val[i][j]=0.0;
		    }
		    else
		    {
			fold.val[i][j]=static_cast<double>(iret);
		    }
		    if(SEISPP_verbose)
                    {
                        if(iret>0)
			    cout << "Grid point ("<<i<<", "<<j
				<<") has stack fold="<<(int)(fold.val[i][j])
				<<endl;
                        else
                        {
                            cout << "Grid point ("<<i<<", "<<j
                                <<") was not processed:  ";;
                            switch(iret)
                            {
                            case (-1):
			    case (0):
                                cout << "stack count below threshold."<<endl;
                                break;
                            case (-2):
                                cout << "centroid of stations outside cutoff"<<endl;
                                break;
                            case (-3):
                                cout << "sum of weights at zero lag is "
                                    << "below threshold (probably 0)"
                                    <<endl;
                                break;
                            default:
                                cout << endl
                                    << "ERROR.  Unknown return code returned by "
                                    << "pwstack_ensemble ="<<iret<<endl
                                    << "This is a coding error that needs to be fixed"
                                    << "Check the source code"<<endl;
                                exit(-1);;
                            }
                        }
                    }
                } catch (SeisppError& serr)
                {
                    serr.log_error();
                    cerr << "Ensemble " << event_number 
                        << " data skipped" << endl;
                    cerr << "Pseudostation grid point indices (i,j)="
                        << "("<<i<<","<<j<<")"<<endl;
                }
            }
	    dfh.save_slowness_vectors(svm,ugrid);
	    char fsbuf[64];
	    sprintf(fsbuf,"%s_%ld",fieldnamebase.c_str(),evid);
	    string fieldname(fsbuf);
            /* This uses the new file based save metod */
            try {
                fold.save(fieldname,fielddir);
            }catch  (exception& e)
            {
                cerr << "Error saving fold data to file "<<fieldname
                    << " in directory "<<fielddir<<endl
                    << "Bundering on as this is not critical data."
                    <<endl;
            }
        }
    } catch (SeisppError& err)
    {
        err.log_error();
    }
    // The GCLgrid library db routines throw simple int exceptions
    // This should only be entered if the GCLgrid constructor fails.
    catch (exception& stdexc)
    {
        cerr << stdexc.what()<<endl;
    }
}
