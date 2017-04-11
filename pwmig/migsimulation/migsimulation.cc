#include <fstream>
#include <vector>
#include <set>
#include "stock.h"
#include "seispp.h"
#include "Metadata.h"
#include "SimplePSPrimarySynthetic.h"
#include "PointSourcePSSynthetic.h"
#include "StaVariableLayeredSynthetic.h"
#include "filter++.h"
#include "SimpleWavelets.h"
#include "VelocityModel_1d.h"
#include "vectorcls.h"

enum SyntheticType {SIMPLE,POINTSOURCE,CONST_VEL_LAYERED};
using namespace std;
using namespace SEISPP;
bool SEISPP::SEISPP_verbose(false);
SyntheticSeismogram *CreateSimpleGenerator(Metadata& control, Pf *pf)
{
    try {
        Tbl *t;
        t=pfget_tbl(pf,const_cast<char *>("interfaces"));
        vector<double> d,a;
        int i;
        for(i=0;i<maxtbl(t);++i)
        {
            char *line;
            line=(char *)gettbl(t,i);
            double depth,amp;
            sscanf(line,"%lf%lf",&depth,&amp);
            d.push_back(depth);
            a.push_back(amp);
        }
        SimplePSPrimarySynthetic *result = new SimplePSPrimarySynthetic(control,d,a);
        return result;
    } catch (...) {throw;};
}
SyntheticSeismogram *CreateStaVariableConstVelGenerator(Metadata& p)
{
    const string base_error("CreateStaVariableConstVelGenerator procedure:  ");
    try{
        double tsigma=p.get_double("tsigma");
        double wlevel=p.get_double("water_level");
        bool useRF=p.get_bool("receiver_function_mode");
        bool single_model_file("use_single_model_file");
        StaVariableLayeredSynthetic *result;
        if(single_model_file)
        {
            string modfile=p.get_string("model_file_name");
            result=new StaVariableLayeredSynthetic(modfile,
                    tsigma,wlevel,useRF);
        }
        else
        {
            string listfile=p.get_string("model_list_file");
            list<string> filelist;
            try {
                ifstream in;
                in.open(listfile.c_str(),ios::in);
                if(in.fail())
                    throw SeisppError(base_error + "Open failed on file "
                            +listfile
                            +"\nThis file should contain a list of velocity model files to read");
                char fname[256];
                while(in.getline(fname,256))
                {
                    filelist.push_back(string(fname));
                }
            } catch(...){throw;};
            result=new StaVariableLayeredSynthetic(filelist,
                    tsigma,wlevel,useRF);
        }
        return result;
    } catch (...) {throw;};
}

SyntheticSeismogram* CreatePointsourceGenerator(Metadata& control,Pf *pf){
//ToDo: implement an interface to initiate an object of Class PointSourcePSSynthetic
//to initiate PointSource obj, a 1D velocity model has to be loaded from the database
//so which velocity model to choose?
int dberror;
try{
Dbptr dbvmodel;
string wavetype;//Svelocity or Pvelocity
string vmodelname;
string fname,txtvmodel;
fname=control.get_string("Velocity_Model_DB");
vmodelname=control.get_string("Velocity_Model");
//wavetype=control.get_string("VModel_paraname");
char *pdbv=new char[fname.size()];
strcpy(pdbv, fname.c_str());
string vmodel_datasource=control.get_string("vmodel_datasource");
txtvmodel=control.get_string("Vmodel_text_file");
if(vmodel_datasource=="database"){
 dberror=dbopen(pdbv, const_cast<char *>("r"), &dbvmodel);
 //dberror=dbopen(&fname[0], "r", &dbvmodel);
 VelocityModel_1d vs1d(dbvmodel, vmodelname, "Svelocity");
 VelocityModel_1d vp1d(dbvmodel, vmodelname, "Pvelocity");
 PointSourcePSSynthetic *result=new PointSourcePSSynthetic(vs1d, vp1d, pf);
 return result;

}
else{
 VelocityModel_1d vs1d(txtvmodel, "plain", "S");
 VelocityModel_1d vp1d(txtvmodel, "plain", "P");
 PointSourcePSSynthetic *result=new PointSourcePSSynthetic(vs1d, vp1d, pf);
 return result;

}

//PointSourcePSSynthetic *result=new PointSourcePSSynthetic(vs1d, vp1d, pf);
//return result;
}
catch(...) {
throw ;
//die(0, "Unable to open the database file. Please verify the file name")
};
}
void build_db_view(DatascopeHandle& dbh, Metadata& control,Pf *pf)
{
    try{
        /* The following is copied exactly from pwstack.  Not very
           elegant, but a way to guarantee the gathers we bring
           in are compatible with pwstack.  If that code changes this
           will need to change too. */
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
		if(SEISPP_verbose) cout << "wfdisc<-arrival join size="
			<<dbh.number_tuples()<<endl;
		dbh.natural_join("assoc");
		dbh.natural_join("origin");
		dbh.natural_join("event");
		dbh.subset("orid==prefor");
		if(SEISPP_verbose) 
                    cout << "wfdisc,<-arrival->assoc->origin->event size="
			<<dbh.number_tuples()<<endl;
                /* Natural join is problematic - use explicity keys
                   so we don't have to run dbfixchanids which requires
                   a sensor table.  */
		//dbh.natural_join("sitechan");
                j1.clear();
		j1.push_back("sta");
                j1.push_back("chan");
		j1.push_back("wfdisc.time::wfdisc.endtime");
                j2.clear();
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
		if(SEISPP_verbose) cout << "Number of event:sta groups="
			<<dbh.number_tuples()<<endl;
	}
	else if(dbviewmode=="use_wfprocess")
	{
		dbh.lookup("event");
		dbh.natural_join("origin");
		dbh.subset("orid==prefor");
		dbh.natural_join("assoc");
		dbh.natural_join("arrival");
		if(SEISPP_verbose) 
		  cout << "Number of rows in event->origin->assoc->arrival view="
		  	<< dbh.number_tuples()<<endl;
		DatascopeHandle ljhandle(dbh);
		ljhandle.lookup("wfprocess");
		ljhandle.natural_join("sclink");
		ljhandle.natural_join("evlink");
		list<string> jk;
		jk.push_back("evid");
		jk.push_back("sta");
		dbh.join(ljhandle,jk,jk);
		if(SEISPP_verbose) 
		  cout << "Number of rows of view after leftjoin of wfprocess view"
		  	<< dbh.number_tuples()<<endl;
		jk.clear();
		jk.push_back("sta");
		dbh.join("site",jk,jk);
		if(SEISPP_verbose) 
		  cout << "Number of rows in final working view="
		  	<< dbh.number_tuples()<<endl;
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
        list<string> group_keys;
        group_keys.push_back("evid");        
        dbh.group(group_keys);        
	if(SEISPP_verbose) 
           cout << "Number of ensembles to process="<<dbh.number_tuples()<<endl;
        dbh.rewind();
    }catch(...) {throw;};
}

set<int> load_eventset(string fname)
{
    ifstream fin;
    fin.open(fname.c_str(),ios::in);
    if(fin.fail()) throw SeisppError(string("load_eventset:  ")
            +"Cannot open file="+fname);
    int val;
    set<int> result;
    if(SEISPP_verbose) cout << "List of Events To Replicate"<<endl;
    while(fin.good())
    {
        fin>>val;
        result.insert(val);
        if(SEISPP_verbose)cout << val<<endl;
    }
    fin.close();
    return(result);
}

void usage()
{
    cerr << "migsimulation dbin dbout [-evf eventlistfile -pf pffile -v]"<<endl;
    exit(-1);
}
SyntheticType parse_syntype(Metadata& type)
{
SyntheticType result;
string strtype;

strtype=type.get_string("Synthetic_method");

    if(strtype=="SIMPLE")
        result=SIMPLE;
    else if(strtype=="POINTSOURCE")
        result=POINTSOURCE;
    else if(strtype=="STAVARIABLE_LAYERED")
        result=CONST_VEL_LAYERED;
    else
        result=SIMPLE;
return result;
}




int main(int argc, char **argv){

    //WARNING WARNING:  frozen for now.  Change when a new synthetic generator is added 
    SyntheticType syntype;//(SIMPLE);
    if(argc<3) usage();
    string dbin_name, dbout_name, pfname("migsimulation");
    dbin_name=string(argv[1]);
    dbout_name=string(argv[2]);
    string evlfile("NONE");
    bool check_evlist(false);
    int i;
    for(i=3;i<argc;++i)
    {
        string strarg(argv[i]);
        if(strarg=="-pf")
        {
            ++i;
            if(i==argc) usage();
            pfname=string(argv[i]);
        }
        else if(strarg=="-evf")
        {
            ++i;
            if(i==argc) usage();
            evlfile=string(argv[i]);
            check_evlist=true;
        }
        else if(strarg=="-v")
            SEISPP_verbose=true;
        else
            usage();
    }
    SyntheticSeismogram *synbase;
    try {
        /* For the present freeze this to be css3.0 schema */
        AttributeMap am("css3.0");
        DatascopeHandle dbhin(dbin_name,true);
        DatascopeHandle dbhout(dbout_name,false);
        DatascopeHandle dbhsc(dbhout);
        dbhsc.lookup("sclink");
        DatascopeHandle dbhev(dbhout);
        dbhev.lookup("evlink");
        dbhout.lookup("wfprocess");
        Pf *pf;
        if(pfread(const_cast<char *>(pfname.c_str()),&pf))
        {
            cerr << "pfread failed for pffile="<<pfname<<endl;
            usage();
        }
        Metadata control(pf);
	syntype=parse_syntype(control);// check the parameter file
        //to decide which simulation method to use
        double tsfull = control.get_double("data_time_window_start");
        double tefull = control.get_double("data_time_window_end");
        TimeWindow data_window(tsfull,tefull);
        /* These should be identical to those defined for pwstack */
        MetadataList station_mdl=pfget_mdlist(pf,"station_metadata");
        MetadataList ensemble_mdl=pfget_mdlist(pf,"ensemble_metadata");
        MetadataList output_mdl=pfget_mdlist(pf,"wfprocess_metadata");
        MetadataList wfdisc_mdl=pfget_mdlist(pf,"wfdisc_metadata");
        /* for now we freeze these names for output */
        vector<string> chanmap;
        chanmap.push_back("E");
        chanmap.push_back("N");
        chanmap.push_back("Z");
        string output_dir=control.get_string("output_directory");
        string dfile_base=control.get_string("output_data_file_base_name");
        /* the next group sets if synthetics are filtered with an
           idealized wavelet or with a standard time invariant filter */
        string wavelet_type=control.get_string("wavelet_type");
        cout << "Using wavelet type = "<<wavelet_type<<endl;
        TimeSeries wavelet;
        double datadt=control.get_double("data_sample_interval");
        // these are not needed for TimeInvariantFilter, but are 
        // needed for gaussian or ricker wavelets 
        int nwsamp = control.get_int("wavelet_length");
        double wavelet_width_parameter
            =control.get_double("wavelet_width_parameter");
        WaveletNormalizationMethod nm=PEAK;
        TimeInvariantFilter *filter;
        if(wavelet_type=="filter")
        {
            /* Use a none definition for filter to do nothing to synthetic
               output */
            string filparam=control.get_string("filter");
            filter=new TimeInvariantFilter(filparam);
        }
        else if(wavelet_type=="gaussian")
        {
            wavelet=gaussian_wavelet(nwsamp,datadt,
                    wavelet_width_parameter,nm);
        }
        else if(wavelet_type=="ricker")
        {
            wavelet=ricker_wavelet(nwsamp,datadt,
                    wavelet_width_parameter,nm);
        }
        else
        {
            cout << "-"<<wavelet_type<<"-"<<endl;
            cout << "wavelet_type parameter="<<wavelet_type
                << " is not supported."<<endl
                << "Must be one of:  filter, gaussian, or ricker"<<endl;
            usage();
        }
        cout << "Wavelet type is set"<<endl;
        /* build event (evid) list if requested */
        set<int> eventset;
        if(check_evlist)
        {
            eventset=load_eventset(evlfile);
        }

        /* Construct the synthetic seismogram generator objects */
        switch(syntype){
            /* For now only one option and it is also default */
         case SIMPLE:
            synbase=CreateSimpleGenerator(control,pf);	
	    break;
         case POINTSOURCE:
	    synbase=CreatePointsourceGenerator(control,pf);
	    break;
         case CONST_VEL_LAYERED:
            synbase=CreateStaVariableConstVelGenerator(control);
            break;
	 default:
            // Necessary because I wanted to not use a pf for this object 
            synbase=CreateSimpleGenerator(control,pf);
        };
        build_db_view(dbhin,control,pf);
        int nevents=dbhin.number_tuples();
        if(SEISPP_verbose) cout << "Processing "<<nevents<<" simulated events"
                            <<endl;
        int rec;
        ThreeComponentSeismogram syndata;
        for(rec=0;rec<nevents;++rec,++dbhin)
        {
            int evid=dbhin.get_int("evid");
            if(check_evlist)
            {
                if(eventset.find(evid)==eventset.end()) continue;
            }
            if(SEISPP_verbose) cout << "Working on evid="<<evid<<endl;
            auto_ptr<ThreeComponentEnsemble> ensemble(
                  new ThreeComponentEnsemble(dynamic_cast<DatabaseHandle&>(dbhin),
                        station_mdl, ensemble_mdl,am));
            ensemble=ArrivalTimeReference(*ensemble,"arrival.time",
                data_window);
            double slat,slon,sz,otime;
            slat=ensemble->get_double("origin.lat");
            slon=ensemble->get_double("origin.lon");
            slat=rad(slat);
            slon=rad(slon);
            sz=ensemble->get_double("origin.depth");
            otime=ensemble->get_double("origin.time");
            //int evid=ensemble->get_int("evid");
            //For now freeze this as taup calculator iasp91
            Hypocenter hypo(slat,slon,sz,otime,string("tttaup"),string("iasp91"));
	    if(syntype==POINTSOURCE)
		((PointSourcePSSynthetic*) synbase)->initPtime4event(hypo);
            vector<ThreeComponentSeismogram>::iterator dptr;
            for(dptr=ensemble->member.begin();
                dptr!=ensemble->member.end();++dptr)
            {
                double rlat,rlon,relev;
		// Ignore data marked data but always log this as a nonfatal erro
		if(!(dptr->live)) 
		{
                    double timestamp=dptr->get_double("arrival.time");
		    cerr << "Deleting data from station="
                        <<dptr->get_string("sta")<<" arrival time tag= "
                        <<strtime(timestamp)<<endl
				<< "Marked dead by constructor.  Likely data problem you need to address"<<endl;
		    continue;
		}
                //DEBUG
                /*cout << dynamic_cast<Metadata &>(*dptr)<<endl;
		cout << dptr->get_string("sta")<<endl;
		cout << "ns="<<dptr->ns<<" dt="<<dptr->dt<<" live="<<dptr->live<<endl;
                */
                rlat=dptr->get_double("site.lat");
                rlon=dptr->get_double("site.lon");
                rlat=rad(rlat);
                rlon=rad(rlon);
                relev=dptr->get_double("site.elev");
                try {
                    syndata = synbase->Compute3C(*dptr,
                            hypo,rlat,rlon,-relev,string(""));
                }catch (SeisppError &serr)
                {
                    cerr << "Error processing station "<<dptr->get_string("sta")
                        << " for event "<< evid<<endl
                        << "No output for this seismogram."<<endl
                        << "Error mesage thrown by synthetic Compute3C method:"
                        <<endl;
                    serr.log_error();
                    continue;
                }
                if(wavelet_type=="filter")
                {
                    filter->apply(syndata);
                }
                else
                {
                    syndata=sparse_convolve(wavelet,syndata);
                }
                /* Data are written to event files in a user 
                   selected directory.  This puts this in the trace
                   headers/metadata area. */
                syndata.put("dir",output_dir);
                stringstream ss;
                ss<<dfile_base<<"_"<<evid;
                syndata.put("dfile",ss.str());
                string sta=dptr->get_string("sta");
                syndata.put("sta",sta);
                /* We have to put syndata back to an absolute time 
                   frame.  Done through arrival.time header field */
                double atime=syndata.get_double("arrival.time");
                syndata.rtoa(atime);
                string algorithm("migsimulation");
                syndata.put("wfprocess.algorithm",algorithm);
                /* This is hack fix for handling data read with wfdisc
                   In that case timetype is not set so we force it always*/
                syndata.put("timetype","a");
                int wfrec=dbsave(syndata,dbhout.db,string("wfprocess"),
                        output_mdl,am);
                int pwfid;
                dbhout.db.record=wfrec;
                pwfid=dbhout.get_int("pwfid");
                /* use datatype for chan here.  May be trouble, but 
                   a starting point */
                string chan=dbhout.get_string("datatype");
                dbhev.append();
                dbhev.put("pwfid",pwfid);
                dbhev.put("evid",evid);
                dbhsc.append();
                dbhsc.put("sta",sta);
                dbhsc.put("chan",chan);
                dbhsc.put("pwfid",pwfid);
                /* now deal with wfdisc.  This assumes
                 dbsave does a lookup on dbhout for wfdisc.
                 Note wfdisc data will be interleaved with wfprocess
                 data all in one large event file*/
                wfrec=dbsave(syndata,dbhout.db,string("wfdisc"),
                        wfdisc_mdl,am,chanmap,true);
            }
        } 
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(dmatrix_error& derr)
    {
        derr.log_error();
    }

}

