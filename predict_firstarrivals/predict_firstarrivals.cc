/*
Compute first arrival P or S and save the arrival time to the arrival table. This code was motivated by the ideal of get_predicted_Parrivals.
Since that Perf script is too slow to process large dataset, I decided to write a similar C++ code to replace it.

Capbilities:
1. compute first arrivals
2. save to db
3. option to read in and save to text files
4. designed for processing large datasets

BRTT Antelope libraries used:
pphasetime()
sphasetime()

*/
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include "tt.h"
#include "db.h"
#include "stock.h"
#include "tt.h"
#include "seispp.h"
#include "SeisppKeywords.h"
#include "dbpp.h"
#include "Metadata.h"

using namespace std;
using namespace SEISPP;

const string MDL_SITE("mdlist_site"), MDL_EVENT("mdlist_event"), 
			MDL_ARRIVAL("mdlist_arrival");
			
void history_current()
{
	cout<<"Feb 11, 2019: created"<<endl
	<<endl;
}

const string csversion("v1.0");

void version()
{
	cerr <<"< version "<<csversion<<" > 2/11/2019"<<endl;
}
void author()
{
	cerr <<endl<<"Xiaotao Yang, University of Massachusetts Amherst"<<endl<<endl;
}
void usage_message()
{
    version();
    cerr << "predict_firstarrivals eventdb sitedb arrivaldbout [-ss subset_condition][-v|V][-h|H]"<<endl;
    cerr << "** Use -h to print out detailed explanations on the options."<<endl;
    author();
}
void help()
{
	usage_message();
	cout<<"TBA"<<endl;
    exit(0);
}
void usage()
{
	usage_message();
    exit(-1);
}

bool SEISPP::SEISPP_verbose(false);

MetadataList generate_mdlist(string mdltype, bool use_arrival_data=false, bool use_netmag_table=false,
			bool use_decon_in_editing=false)
{
	try
	{
		MetadataList mdlist;
		Metadata_typedef metadata;
	
		if(mdltype=="ensemble")
		{
			metadata.tag="sta";
			metadata.mdt=MDstring;
			mdlist.push_back(metadata);
		}
		else if(mdltype=="wfprocessin")
		{
			metadata.tag="pwfid"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="dir"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="dfile"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="time"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="endtime"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="nsamp"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="samprate"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="datatype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="timetype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="foff"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="wfprocess.algorithm"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="evid"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="sta"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="chan"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			if(use_decon_in_editing)
			{
				metadata.tag="decon.nspike"; metadata.mdt=MDint; mdlist.push_back(metadata);
				metadata.tag="decon.rawsnr"; metadata.mdt=MDreal; mdlist.push_back(metadata);
				metadata.tag="decon.averamp"; metadata.mdt=MDreal; mdlist.push_back(metadata);
				metadata.tag="decon.epsilon"; metadata.mdt=MDreal; mdlist.push_back(metadata);
				metadata.tag="decon.niteration"; metadata.mdt=MDint; mdlist.push_back(metadata);
				metadata.tag="decon.peakamp"; metadata.mdt=MDreal; mdlist.push_back(metadata);
				metadata.tag="decon.chan"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			}
			if(use_arrival_data)
			{
				metadata.tag="atime"; metadata.mdt=MDreal; mdlist.push_back(metadata);
				metadata.tag="arrival.sta"; metadata.mdt=MDstring; mdlist.push_back(metadata);
				metadata.tag="assoc.esaz"; metadata.mdt=MDreal; mdlist.push_back(metadata);
				metadata.tag="assoc.seaz"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			}
			if(use_netmag_table)
			{
				metadata.tag="magtype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
				metadata.tag="magnitude"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			}
		}
		else if(mdltype=="wfprocessout")
		{
			metadata.tag="pwfid"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="time"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="endtime"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="dir"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="dfile"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="datatype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="timetype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="samprate"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="nsamp"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="foff"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="wfprocess.algorithm"; metadata.mdt=MDstring; mdlist.push_back(metadata);
		}
		else if(mdltype=="wfdiscin")
		{
			metadata.tag="wfdisc.chan"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="wfdisc.time"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="wfdisc.wfid"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="wfdisc.chanid"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="wfdisc.jdate"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="wfdisc.endtime"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="wfdisc.nsamp"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="wfdisc.samprate"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="wfdisc.calib"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="wfdisc.calper"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="wfdisc.instype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="wfdisc.segtype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="wfdisc.datatype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="wfdisc.clip"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="wfdisc.dir"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="wfdisc.dfile"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="wfdisc.foff"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="wfdisc.commid"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="sta"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="chan"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="time"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="nsamp"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="samprate"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			if(use_arrival_data)
			{
				metadata.tag="atime"; metadata.mdt=MDreal; mdlist.push_back(metadata);
				metadata.tag="arrival.sta"; metadata.mdt=MDstring; mdlist.push_back(metadata);
				metadata.tag="assoc.esaz"; metadata.mdt=MDreal; mdlist.push_back(metadata);
				metadata.tag="assoc.seaz"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			}
			if(use_netmag_table)
			{
				metadata.tag="magtype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
				metadata.tag="magnitude"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			}
		}
		else if(mdltype=="wfdiscout")
		{
			metadata.tag="sta"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="chan"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="time"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="endtime"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="nsamp"; metadata.mdt=MDint; mdlist.push_back(metadata);
			metadata.tag="samprate"; metadata.mdt=MDreal; mdlist.push_back(metadata);
			metadata.tag="datatype"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="dir"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="dfile"; metadata.mdt=MDstring; mdlist.push_back(metadata);
			metadata.tag="foff"; metadata.mdt=MDint; mdlist.push_back(metadata);
		}
		else
		{
			cerr<<"ERROR in generate_mdlist(): wrong metadata type!"<<endl;
			exit(-1);
		}
	
		return(mdlist);
	}catch(...) {throw;};
}
//define the priority list in matching channels
list<string> get_chanlist()
{
	try
	{
		list<string> chanlist;
		chanlist.push_back("BHZ");
		chanlist.push_back("HHZ");
		chanlist.push_back("EHZ");
		return(chanlist);
	}catch(...) {throw;};
}
/*==================================================================
//==================================================================
//====================== Main program ==============================
//==================================================================
//==================================================================
*/
int main(int argc, char **argv)
{
    const string errorbase("Error in predict_firstarrivals: ");
    switch(argc)
    {
    	case 1:
    		usage();
    		break;
    	case 2:
    		string sarg(argv[1]);
    		if(sarg=="-h" || sarg=="-H")
    			help();
    		else if(sarg=="-history")
        	{	history_current();
        		exit(0);}
        	else
    			usage();
    		break;
    }
    string eventdb_name;
	string sitedb_name;
	string outdb_name;
	if (argc==3)
		usage();
	else if (argc==4)
	{
		eventdb_name=argv[1];
		sitedb_name=argv[2];
		outdb_name=argv[3];
	}
    bool save_to_text_file(false);
    string outfname("predicted_arrivals_summary.txt");

    bool apply_site_subset(false), apply_origin_subset(false);
    bool set_continue_mode_by_default(false);
    bool turn_on_continue_mode(true);

    	//unless review-mode and GUIoff are on, post-edit-FA will be computed.
    string site_subset_condition(""),origin_subset_condition("");
	/* quotes needed around subsets because sometimes
	station names start with numbers that confuse datascope */
	const string quote("\"");
    int i,j;
    for(i=4;i<argc;++i)
    {
        string sarg(argv[i]);
        
        if(sarg=="-continue")
        {
        	set_continue_mode_by_default=true; 
        		/*
        		//if this is true, the program will continue working on next station 
        		//when output table is not empty. otherwise, it will ask the user to 
        		//choose either continue on the next station or quit.
        		// Xiaotao Yang
        		*/
        }
        else if(sarg=="-ss") //site table subset
        {
            ++i;
            if(i>=argc) usage();
            apply_site_subset=true;
            if(site_subset_condition.length()> 0)
            	site_subset_condition=site_subset_condition+string("&&")+string(argv[i]);
            else
	            site_subset_condition=string(argv[i]);
        }
        else if(sarg=="-os") //origin table subset
        {
            ++i;
            if(i>=argc) usage();
            apply_origin_subset=true;
            if(origin_subset_condition.length()> 0)
            	origin_subset_condition=origin_subset_condition+string("&&")+string("orid > ")
            		+quote+string(argv[i])+quote;
            else
            	origin_subset_condition=string("sta > ")+quote
                	+string(argv[i])+quote;
        }
        else if(sarg=="-v" || sarg=="-V")
        	SEISPP_verbose=true;
        else if(sarg=="-h" || sarg=="-H")
        	help();
        else if(sarg=="-history")
        {	history_current();
        	exit(0);}
        else
            usage();
    }
    
    double  distance=0.0, depth=0.0, phtimes=0.0;
	
	/****************************************************************
		*****************************************************************
		*******<<<<<<<<<<<< READ IN CONTROL PARAMETERS >>>>>>>>>>********
		*****************************************************************
		*****************************************************************
		*/
	string zchan="BHZ"; // make sure this is consistent with the real channel for each site-chan pair
	

		/****************************************************************
		*****************************************************************
		*********<<<<<<<<<<<< CHECK OUTPUT DATABASE >>>>>>>>>>***********
		*****************************************************************
		*****************************************************************
		*/

	DatascopeHandle dbout(outdb_name,false);
	DatascopeHandle dbharrival(dbout);

	/* First check that the output wfdisc is empty */

	dbout.lookup("arrival");
	cout << "Writing results to [ "<<outdb_name<<" ]"<<endl;
	if (dbout.number_tuples()>0)
	{	cout <<"Number of existing rows in the output database table is "
		<<dbout.number_tuples()<<endl;
	}	
	/****************************************************************
	*****************************************************************
	*******<<<<<<<<<<< BUILD WORKING VIEWS >>>>>>>>>********
	*****************************************************************
	*****************************************************************
	*/
	AttributeMap am("css3.0");
	/* Open the in and out database */
	DatascopeHandle dbin_event(eventdb_name,true);
	DatascopeHandle dbin_site(sitedb_name,true);
	
	dbin_event.lookup("origin");

	if (dbin_event.number_tuples()>0)
	{
		cout <<"Number of rows in the origin table before subset is "
			<<dbin_event.number_tuples()<<endl;
	}
	else
	{
		cerr << errorbase << "No events found in the origin table of database [ "<<eventdb_name<<" ] 1"<<endl;
        exit(-1);
	}
	if (apply_origin_subset)
	{
		dbin_event.subset(string(origin_subset_condition));
		cout << "Applying subset for origin table: " << origin_subset_condition << endl;
		cout <<"Number of rows in the origin table after subset is "
			<<dbin_event.number_tuples()<<endl;
	}
	
	long nevent=0;
	nevent=dbin_event.number_tuples();
	cout <<"Number of origins to processes: "<<nevent<<endl;
	
	
	dbin_site.lookup("sitechan"); //use sitechan table here to consider the channel updates at different time periods
	if (apply_site_subset)
	{
		dbin_site.subset(string(site_subset_condition));
		cout << "Applying subset for origin table: " << site_subset_condition << endl;
	}
	

	if (dbin_site.number_tuples()>0)
	{
		cout <<"Number of rows in the sitechan table is "
			<<dbin_site.number_tuples()<<endl;
	}
	else
	{
		cerr << errorbase << "No sites found in the sitechan table of database [ "<<sitedb_name<<" ] 1"<<endl;
        exit(-1);
	}
	
	list<string> sortkeys;
	sortkeys.push_back("orid");
    dbin_event.sort(sortkeys); // sort the origin table by orid
    double otime, sctime_on, sctime_off; //otime: origin time; sctime_on: site chan ontime; sctime_off: site chan offtime.
    cout <<">>++++++++++++++++++++++++++++++"<<endl;
    
//     for(i=0,dbin_event.rewind();i<nevent;++i,++dbin_event)
	for(i=0,dbin_event.rewind();i<3;++i,++dbin_event)
	{   
		cout <<"Working on orid: "<<dbin_event.get_int("orid")<<" ... "<<i+1<<" of "<<nevent<<""<<endl;
		otime=dbin_event.get_double("time");
		cout << otime << endl;
		
		DatascopeHandle dbin_site_temp(dbin_site);
		dbin_site_temp.lookup("sitechan");
		cout <<"Number of rows in the sitechan table master is: "<<dbin_site_temp.number_tuples()<<endl;
		
		string site_subset_local=string("ondate <= ")+epoch2str(otime,"%Y%j")+string("&& offdate >= ")+epoch2str(otime,"%Y%j");
		
		cout <<"subset site with: "<<site_subset_local<<endl;
		dbin_site_temp.subset(site_subset_local);
		list<string> sortkeys2, groupkeys2;
        sortkeys2.push_back("sta");
        sortkeys2.push_back("chan");
        groupkeys2.push_back("sta");
        dbin_site_temp.sort(sortkeys2);
        DatascopeHandle dbin_site_bkp(dbin_site_temp);
        dbin_site_temp.group(groupkeys2);
		int nsta(dbin_site_temp.number_tuples()),nline(dbin_site_bkp.number_tuples());
		
		cout <<"Number of sites to process: "<<nsta<<endl;
		
		for (j=0,dbin_site_bkp.rewind();j<nline;++j,++dbin_site_bkp)
		{
			string chan;
			chan=dbin_site_bkp.get_string("chan");
			if (chan=="BHZ")
				cout<<chan<<endl;
		}
				
    }
    
	distance=50;
	depth=10;
	phtimes = pphasetime(distance,depth);
        
	printf("%10.3f  %10.3f  %10.3f\n",phtimes,distance,depth);


    return 0;
}










