/*
Compute first arrival P or S and save the arrival time to the arrival table. This code was motivated by the ideal of get_predicted_Parrivals.
Since that Perf script is too slow to process large dataset, I decided to write a similar C++ code to replace it.

Capbilities:
1. compute first arrivals
2. save to db
3. option to read in and save to text files
4. designed for processing large datasets

BRTT Antelope functions used:
    pphasetime()
    sphasetime()
    dist()
    epoch2str()
    yearday()

 Xiaotao Yang, UMass Amherst, Feb 11-12, 2019
*/
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include "tt.h"
#include "db.h"
#include "coords.h"  //dist() function
#include "stock.h"
#include "tt.h"
#include "seispp.h"
#include "dbpp.h"

using namespace std;
using namespace SEISPP;

const double PI(3.141592653);
			
void history_current()
{
	cout<<"Feb 11-12, 2019: created"<<endl
	<<endl;
}

const string csversion("v1.0");

void version()
{
	cerr <<"< version "<<csversion<<" > 2/12/2019"<<endl;
}
void author()
{
	cerr <<endl<<"Xiaotao Yang, University of Massachusetts Amherst"<<endl<<endl;
}
void usage_message()
{
    version();
    cerr << "predict_firstarrivals eventdb sitedb arrivaldbout "<<endl
	 << "      [-ph phasename][-ss site_subset][-os origin_subset][-continue][-v|V][-h|H]"<<endl;
    cerr << "** Use -h to print out detailed explanations on the options."<<endl
         << "** Only mark arrivals on *Z channels."<<endl;
    author();
}
void help()
{
	usage_message();
    cout<<endl
    <<"<<Required: eventdb sitedb arrivaldbout >>"<<endl
        <<endl
        <<"<<Options>>"
	<<"-ph phasename: "<<endl
	<<"                 Specify phase types: P or S. Default is P."<<endl
        <<"-ss site_subset: "<<endl
        <<"                 Subset expression for sitechan table. "<<endl
        <<"                 Subset for vertical channels is built-in in the code."<<endl
        <<"-os origin_subset: "<<endl
        <<"                 Subset condition for origin table."<<endl
        <<"-continue: "<<endl
        <<"                 When arrivaldbout is NOT empty. They program will "<<endl
        <<"                 contrinue adding to it without asking the user! Be cautious!!!"<<endl
        <<"-v|V: "<<endl
        <<"                 Verbose mode."<<endl
        <<"-h|H: "<<endl
        <<"                 This help message."<<endl;
    exit(0);
}
void usage()
{
	usage_message();
    exit(-1);
}

bool SEISPP::SEISPP_verbose(false);

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
	else if (argc>=4)
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
    string phasename("P");
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
        else if(sarg=="-ph") //specify phases: P or S
        {
            ++i;
            if(i>=argc) usage();
            phasename=string(argv[i]);
            if (phasename != "P" && phasename != "S" && phasename != "p" && phasename != "s")
            {
                cerr<<errorbase<<"wrong phase name. Can ONLY be P or S."<<endl;
                exit(0);
            }
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
            	origin_subset_condition=origin_subset_condition+string("&&")+string(argv[i]);
            else
            	origin_subset_condition=string(argv[i]);
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

	string zchan="BHZ";
    long arid(1); //count for arrival id.
    /* ******************************************************************
		*********<<<<<<<<<<<< CHECK OUTPUT DATABASE >>>>>>>>>>***********
		*****************************************************************
		*/
	DatascopeHandle dbout(outdb_name,false);

	/* First check that the output wfdisc is empty */
    cout << "Writing results to [ "<<outdb_name<<" ]"<<endl;
	dbout.lookup("arrival");
	if (dbout.number_tuples()>0)
	{	cout <<"Number of existing rows in the output database table is "
            <<dbout.number_tuples()<<endl;
        if (!set_continue_mode_by_default)
        {
            cerr<<"Output db is NOT empty. Specify -continue OR clean up the outdb manually. Quit now."<<endl;
            exit(-1);
        }
        else
        {
            if (SEISPP_verbose)
                cout <<"!!! Warning: Continue mode is on! Adding to the table WITHOUT checking for duplicates!!!"<<endl;
            list<string> sortkeysout;
            sortkeysout.push_back("arid");
            dbout.sort(sortkeysout);
            dbout.set_record(dbout.number_tuples()-1);
//            cout<<dbout.current_record()<<endl;
            arid=dbout.get_int("arid");
            if (SEISPP_verbose) cout<<"last arid: "<<arid<<endl;
            ++arid;
        }
	}
    try
    {
        /****************************************************************
        *******<<<<<<<<<<< BUILD WORKING VIEWS >>>>>>>>>********
        *****************************************************************
        */
//        AttributeMap am("css3.0");
        /* Open the input database */
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
        
        dbin_site.lookup("sitechan");
        if (dbin_site.number_tuples()>0)
        {
            cout <<"Number of rows in the original sitechan table is: "
            <<dbin_site.number_tuples()<<endl;
        }
        //use sitechan table here to consider the channel updates at different time periods
        //first subset to include only vertical channels for each station.
        if (apply_site_subset)
            site_subset_condition=site_subset_condition+string(" && ")+string("chan=~/.*Z/");
        else
            site_subset_condition=string("chan=~/.*Z/");
        
        dbin_site.subset(string(site_subset_condition));
        cout << "Applying subset for sitechan table: " << site_subset_condition << endl;

        if (dbin_site.number_tuples()>0)
        {
            cout <<"Number of rows in the sitechan table after subset is "
            <<dbin_site.number_tuples()<<endl;
        }
        else
        {
            cerr << errorbase << "No sites found in the sitechan table of database [ "
                <<sitedb_name<<" ] "<<endl;
            exit(-1);
        }

        list<string> sortkeys,sortkeys_site;
        sortkeys.push_back("orid");
        dbin_event.sort(sortkeys); // sort the origin table by orid
        sortkeys_site.push_back("sta"); //sort the sitechan table by sta and then by chan
        sortkeys_site.push_back("chan");
        dbin_site.sort(sortkeys_site);
        
        double otime, sctime_on, sctime_off; //otime: origin time; sctime_on: site chan ontime; sctime_off: site chan offtime.
        double olat,olon,odepth,slat,slon; //locations of the event and the site.
        double  distance(0.0), phtimes(0.0),ftemp(0.0);
//        long newrowidx;
        cout <<">>++++++++++++++++++++++++++++++"<<endl;
        
    //     for(i=0,dbin_event.rewind();i<nevent;++i,++dbin_event)
        for(i=0,dbin_event.rewind();i<nevent;++i,++dbin_event)
        {
            cout <<"> Working on orid: "<<dbin_event.get_int("orid")<<" ... "<<i+1<<" of "<<nevent<<""<<endl;
            otime=dbin_event.get_double("time");
            olat=dbin_event.get_double("lat");
            olon=dbin_event.get_double("lon");
            odepth=dbin_event.get_double("depth");
    //        cout << strtime(otime) << endl;
            
//            DatascopeHandle dbin_site_temp(dbin_site);
            dbin_site.lookup("sitechan");
            dbin_site.natural_join("sitechan","site");
            
            if (SEISPP_verbose) {
                dbin_site.subset(site_subset_condition);
                cout <<"Number of sites in the sitechan table master is: "
                    <<dbin_site.number_tuples()<<endl;
            }
            //only compute for sites operating during the period covering the event.
            string site_subset_local=site_subset_condition+string(" && ")+string("ondate <= ")+epoch2str(otime,"%Y%j")+string("&& offdate >= ")+epoch2str(otime,"%Y%j");
            
            if (SEISPP_verbose) cout <<"subset site with: "<<site_subset_local<<endl;
            dbin_site.subset(site_subset_local);
            
            int nsta(dbin_site.number_tuples());
            
            cout <<"Number of sites to process (may not be unique): "<<nsta<<endl;
            
            if (nsta>0)
            {
                for (j=0,dbin_site.rewind();j<nsta;++j,++dbin_site)
                {
                    zchan=dbin_site.get_string("chan");
    //                if (SEISPP_verbose) cout<<zchan<<endl;
                    slat=dbin_site.get_double("lat");
                    slon=dbin_site.get_double("lon");
                    
                    //get great circle distance.
                    dist(olat,olon,slat,slon,&distance,&ftemp);
                    //convert to degrees from radians.
                    distance=distance*180.0/PI;
                    
                    //get phase time.
                    if (phasename == "P" || phasename == "p")
                    {
                        phasename="P";
                        phtimes = pphasetime(distance,odepth);
                    }
                    else if(phasename == "S" || phasename == "s")
                    {
                        phasename="S";
                        phtimes = sphasetime(distance,odepth);
                    }
                    
        //            printf("%10.3f  %10.3f  %10.3f\n",phtimes,distance,odepth);
        //            cout<<phtimes+otime<<endl;
        //            cout<<epoch2str(phtimes+otime,"%Y%j")<<endl;
                    dbout.lookup("arrival");
                    dbout.append();
//                    cout<<newrowidx<<endl;
                    dbout.put("sta",dbin_site.get_string("sta"));
                    dbout.put("chan",zchan);
                    dbout.put("time",phtimes+otime);
                    dbout.put("arid",arid);
                    dbout.put("jdate",yearday(phtimes+otime));
                    dbout.put("iphase",phasename);
                    dbout.put("auth","predicted");
                    
                    ++arid;
                }
            }
        }
        
        dbout.close();
        dbin_site.close();
        dbin_event.close();
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(exception& stdexcept)
    {
        cerr << "Exception thrown:  "<<stdexcept.what()<<endl;
    }
}
