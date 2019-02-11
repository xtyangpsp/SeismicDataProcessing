#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include "ray1d.h"  //libseispp library
#include "tt.h"
#include "seispp.h"
#include "SeisppKeywords.h"
#include "dbpp.h"
#include "Metadata.h"
using namespace std;
using namespace SEISPP;
//const string evidkey("eventid");   // not evid to avoid collision
const string tracetype_3c("3C"), tracetype_1c("1C");
const string MDL_ENSEMBLE("mdlist_ensemble"), MDL_WFDISCIN("mdlist_wfdisc_in"), 
			MDL_WFDISCOUT("mdlist_wfdisc_out"), MDL_WFPROCESSIN("mdlist_wfprocess_in"),
			MDL_WFPROCESSOUT("mdlist_wfprocess_out");
			
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
    cerr << "RFeditor eventdb sitedb arrivaldbout [-ss subset_condition][-v|V][-h|H]"<<endl;
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

/* This function is used by both the TimeSeries and Three_Component
versions of dbsave below.  It builds a database row from a metadata
object, which is produced in both cases by  casting up to Metadata,
and pushing attributes out driven by the list, mdl, and the 
namespace map, am.  

Arguments:
	md = Metadata object containing attributes to be written to 
		database
	db = Antelope Dbptr.  It MUST point at a valid single row that
		is to contain the attributes to be copied there.  
	table = name of table these attributes are being written to
		(needed for consistency check).
	mdl = defines names to be extracted from metadata to
		write to database output
	am = AttributeMap object defining internal to external namespace
		mapping

The basic algorithm is:
	for each element of mdl
		if mdl-> am -> dbtable_name == table
			save
		else
			skip
	end foreach
*/
//Copied from readwrite.cc of SEISPP libs by Xiaotao Yang.
void save_metadata_for_object2(Metadata& md,
	Dbptr db,
		string table,
			MetadataList& mdl, 
				AttributeMap& am)
		throw(SeisppError)
{
	MetadataList::iterator mdli;
	map<string,AttributeProperties>::iterator ami,amie=am.attributes.end();
	map<string,AttributeProperties> aliasmap;
	string cval;
	const string base_message("dbsave->save_metadata_for_object2:  ");
//DEBUG
	for(mdli=mdl.begin();mdli!=mdl.end();++mdli)
	{
		double dval;
		long ival;
		string mdkey;
		if(am.is_alias((*mdli).tag))
		{
			try {
				aliasmap=am.aliases((*mdli).tag);
			} catch(SeisppError& serr){throw serr;};
			ami=aliasmap.find(table);
			if(ami==aliasmap.end())
			{
				dbmark(db);
				throw SeisppError(base_message
				 + string("Alias name=")
				 + (*mdli).tag
				 + string(" is not associated with table=")
				 + table
				 + string("\nVerify output specification against schema") );
			}
			mdkey=(*mdli).tag;   // in this case the alias is the key 
		}
		else
		{
			mdkey=(*mdli).tag;
			ami = am.attributes.find(mdkey);
			if(ami==amie) 
			{
				dbmark(db);
				throw SeisppError(
					string("Required attribute ")
					+(*mdli).tag
					+string(" cannot be mapped to output namespace"));
			}
			if( (ami->second.db_table_name) != table)
			{
				dbmark(db);
				throw SeisppError( 
					string("dbsave (database table mismatch): attribute ")
					+ ami->second.db_attribute_name
					+ string(" is tagged with table name ")
					+ ami->second.db_table_name
					+ string("expected to find ")
					+ table);
			}
			/* In this case the key we use the name from the Attribute map as the key */
			mdkey=ami->second.internal_name;
		}
		try {
			switch(ami->second.mdt)
			{
			case MDint:
				if(ami->second.is_key)
				{
					ival = dbnextid(db,
					  const_cast<char *>
					   (ami->second.db_attribute_name.c_str()) );
					if(ival<0)throw SeisppError(
					  	string("dbsave:  ")
						+ ami->second.db_attribute_name
						+ string(" is defined as integer key for table ")
						+ ami->second.db_table_name
						+ string(" but dbnextid failed") );
					
				}
				else
					ival = md.get_long(mdkey);
				dbputv(db,0,ami->second.db_attribute_name.c_str(),
					ival,NULL);
				// In this case we need to push this back to metadata
				// so it can be used downstream
				md.put(ami->second.db_attribute_name,ival);
				break;
			case MDreal:
				dval = md.get_double(mdkey);
				dbputv(db,0,ami->second.db_attribute_name.c_str(),
					dval,NULL);
				break;
			case MDstring:
				cval = md.get_string(mdkey);
				dbputv(db,0,ami->second.db_attribute_name.c_str(),
					cval.c_str(),NULL);
				break;
			case MDboolean:
				// treat booleans as ints for external representation
				// This isn't really necessary as Antelope
				// doesn't support boolean attributes
				if(md.get_bool(mdkey))
					ival = 1;
				else
					ival = 0;
				dbputv(db,0,ami->second.db_attribute_name.c_str(),
					ival,NULL);
				break;
				
			case MDinvalid:
				cerr << "dbsave: database attribute "
					<< ami->second.db_attribute_name
					<< " was marked as invalid\n"
					<< "Data for this attribute not saved"
					<< endl;
				break;
			
			default:
				cerr << "dbsave: database attribute "
					<< ami->second.db_attribute_name
					<< " has unknown data type\n"
					<< "Data for this attribute not saved"
					<< endl;
			}
	
		}
		catch (MetadataGetError& mderr)
		{
			mderr.log_error();
			throw SeisppError(
			    string("dbsave object failure from problem in metadata components"));
		}
	}
}

bool check_continue_mode(bool set_continue_mode_by_default,string laststation)
{
	try
	{
		bool turn_on_continue_mode(false);
		string ques;
		if(set_continue_mode_by_default) 
		{
			turn_on_continue_mode=true;
			return(turn_on_continue_mode);
		}
		else
		{
			cout<<"!!! Output table is not empty. Last station you worked on = [ "
				<<laststation<<" ]."<<endl
				<<"> Continue working on next station? (y/r/n) "<<endl
				<<"	y: continue;"<<endl
				<<"	r: re-do from the beginning. Caution for duplicate rows! "<<endl
				<<"		Suggest choose this only under review mode;"<<endl
				<<"	n: quit the program."<<endl<<"> Your choice: ";
			cin>>ques;
			if(ques == "y" || ques=="Y") turn_on_continue_mode=true;
			else if(ques == "r" || ques == "R") 
			{
				turn_on_continue_mode=false;
				cout<<"!!! CAUTION: if you are not under review mode, "<<endl
					<<"output wfdisc table will have duplicate rows!"<<endl;
				return(turn_on_continue_mode);
			}
			else
			{
				cout<<"Quited. "<<endl<<"Please clean up the output db "
					<<"and the data directory."<<endl;
				exit(-1);
			}					
		}
		return(turn_on_continue_mode);
	}catch(...){throw;};
}
//mdltype: ensemble, wfprocessin, wfdiscin, wfprocessout, wfdiscout
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

/*==================================================================
//==================================================================
//====================== Main program ==============================
//==================================================================
//==================================================================
*/
int main(int argc, char **argv)
{
    const string errorbase("Error in RFeditor: ");
    const string logfilename("RFeditor.log");
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
        	else if(sarg=="-history2")
        	{	history_old();
        		history_current();
        		exit(0);}
    		else 
    			usage();
    		break;
    }
    
    string dbin_name(argv[1]);
    string dbout_name(argv[2]);
    string outdir("RFDataEdited");
    string pfname("RFeditor");
    string FA_filename("-");
    bool save_edit_summary_to_file(false);
    string editsummaryfname("TraceEditSummary.txt");
    //string stacktype("r");  // set default stacktype as RobustSNR. Xiaotao Yang
    string ques;
    bool apply_subset(false);
    bool review_mode(false);
    bool set_continue_mode_by_default(false);
    bool turn_on_continue_mode(true);
    bool datatype3c(false);
    bool GUIoff(false);
    //program will write FirstArrival information into 
    bool get_FA(false);
    bool pre_edit_FA(false),post_edit_FA(false); 
    	//unless review-mode and GUIoff are on, post-edit-FA will be computed.
    string subset_condition("");
    string filterspec;
	/* quotes needed around subsets because sometimes
	station names start with numbers that confuse datascope */
	const string quote("\"");
    int i;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        // read in for stack type
        if(sarg=="-rm" || sarg=="--review-mode")
        {
        	review_mode=true;
        	/*remind user of the review_mode. Xiaotao Yang*/
			if(review_mode)
			{
				cout<<"Warning: review mode is turned on! All edits will be dropped!"<<endl
					<<"> Continue? (y/n) ";
				cin>>ques;
				if(ques != "y" && ques!="Y")
				{
					cout<<"Quited."<<endl;
					exit(-1);
				}
			}
		}
        else if(sarg=="--gui-off" || sarg=="-go")
        	GUIoff=true;
        else if(sarg=="--first-arrival " || sarg=="-fa")
        {
        	get_FA=true,pre_edit_FA=false,post_edit_FA=true;
        	++i;
            if(i>=argc) usage();
            FA_filename=string(argv[i]);
        }
        else if(sarg=="-continue")
        {
        	set_continue_mode_by_default=true; 
        		/*
        		//if this is true, the program will continue working on next station 
        		//when output table is not empty. otherwise, it will ask the user to 
        		//choose either continue on the next station or quit.
        		// Xiaotao Yang
        		*/
        }
        else if(sarg=="-pf")
        {
            ++i;
            if(i>=argc) usage();
            pfname=string(argv[i]);
        }
        else if(sarg=="-d")
        {
        	++i;
            if(i>=argc) usage();
            outdir=string(argv[i]);
        }
        else if(sarg=="-tredit")
        {
        	++i;
            if(i>=argc) usage();
            save_edit_summary_to_file=true;
            editsummaryfname=string(argv[i]);
        }
        else if(sarg=="-ss")
        {
            ++i;
            if(i>=argc) usage();
            apply_subset=true;
            if(subset_condition.length()> 0)
            	subset_condition=subset_condition+string("&&")+string(argv[i]);
            else
	            subset_condition=string(argv[i]);
        }
        else if(sarg=="-laststa")
        {
            ++i;
            if(i>=argc) usage();
            apply_subset=true;
            if(subset_condition.length()> 0)
            	subset_condition=subset_condition+string("&&")+string("sta > ")
            		+quote+string(argv[i])+quote;
            else
            	subset_condition=string("sta > ")+quote
                	+string(argv[i])+quote;
        }
        else if(sarg=="-v" || sarg=="-V")
        	SEISPP_verbose=true;
        else if(sarg=="-h" || sarg=="-H")
        	help();
        else if(sarg=="-history")
        {	history_current();
        	exit(0);}
        else if(sarg=="-history2")
        {	history_old();
        	history_current();
        	exit(0);}
        else
            usage();
    }
    
    //Review Mode must be off in GUIoff mode.
    if(review_mode && GUIoff && !get_FA)
    {
    	cerr<<errorbase<<"Can't turn on Review-Mode while under GUIoff mode "
    		<<"unless [-fa] option is used."<<endl;
    	usage();
    }
    else if(review_mode && GUIoff && get_FA && SEISPP_verbose)
    {
    	cout<<"!!!Caution: review-mode, GUIoff-mode and first-arrival mode are all turned on."<<endl
    		<<"This will result in ONLY writing out first arrival information without doing editing!"<<endl;
    	pre_edit_FA=true, post_edit_FA=false;
    }
    /*
    if(post_edit_FA) 
    	cout<<"!!!Caution: post-edit FA will ONLY be saved when "
    		<<"editing on radial component! Ignored otherwise!"<<endl;
    */
    if(SEISPP_verbose && get_FA)
    	cout<<"Get-FirstArrival mode is on. FA information will be written into [ "<<FA_filename<<" ]."<<endl;
    /*Standard way to open Antelope pf file.*/
    Pf *pf;
    if(pfread(const_cast<char *>(pfname.c_str()),&pf))
    {
        cerr << "pfread failed on pf file named "<<pfname<<endl;
        exit(-1);
    }
    try {
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

		DatascopeHandle dbout(dbout_name,false);
		DatascopeHandle dbharrival(dbout);

		/* First check that the output wfdisc is empty */

		dbout.lookup("arrival");
		cout << "Writing results to "<<dbout_name<<endl
			<<"Number of existing rows in the output database table is "
			<<dbout.number_tuples()<<endl;

		/****************************************************************
		*****************************************************************
		*******<<<<<<<<<<< BUILD WORKING VIEWS >>>>>>>>>********
		*****************************************************************
		*****************************************************************
		*/
        AttributeMap am("css3.0");
        /* Open the in and out database */
        DatascopeHandle dbinevent(dbin_name,true);
        if(use_wfdisc_in) 
        {
        	//dbin.lookup("wfdisc"); 
        		//read in arrival for wfdisc need to be added later.
        	if(use_arrival_data)
        	{
				dbin.lookup("event");
				dbin.natural_join("origin");
				dbin.subset("orid==prefor");
				dbin.natural_join("assoc");
				dbin.natural_join("arrival");
				string phase_subset;
				phase_subset="phase=~/P/";
				dbin.subset(phase_subset);
				if(SEISPP_verbose) 
					cout<< "Size of catalog view ="<< dbin.number_tuples()<<endl;
				list<string> j1,j2;
				j1.push_back("sta");
				j1.push_back("wfdisc.time::wfdisc.endtime");
				j2.push_back("sta");
				j2.push_back("arrival.time");
				dbin.leftjoin(string("wfdisc"),j1,j2);
			
				if(SEISPP_verbose) 
					cout<< "Size of working view after joining with wfdisc ="
						<< dbin.number_tuples()<<endl;
			}
			else
			{
				dbin.lookup("wfdisc");
				if(SEISPP_verbose) 
					cout<< "Size of working view in wfdisc ="
						<< dbin.number_tuples()<<endl;
			}
        }
        else
        //evlink, sclink tables are required when using wfprocess table as input.
        //on the other hand, decon table is optional in case the user doesn't have decon table.
        {
        	if(SEISPP_verbose) cout<<"Building waveform view ..."<<endl;
        	dbin.lookup("wfprocess");
        	dbin.natural_join("evlink");
        	dbin.natural_join("sclink");
        	dbin.db.record=0;
        	string datatype=dbin.get_string("datatype");
        	if(datatype=="c3")	datatype3c=true;
        	else
        		save_3C_data=false; //overwrite the value read from pf. 
        	if(use_decon_in_editing)
        	{	
        		dbin.natural_join("decon");
				/*
				//3c datatype uses one signle row to store 3 component data.
				//not the correct way but reasonable for trace editing (since we only edit radial
				//traces) here we set the subset condition to subset the view after join with decon
				//resulting in only R chan. decon attributes for R chan will be read in as the 
				//metadata for all three components.
				*/
				
				if(datatype3c)
				{
					string sstring=string("decon.chan==")+quote+string(edit_chan_code)+quote;
					dbin.subset(string(sstring));
					/*
					if(subset_condition.length()> 0)
						subset_condition=subset_condition+string("&&")+string("decon.chan==")
							+quote+string(rchan)+quote;
					else
						subset_condition=string("decon.chan==")+quote+string(rchan)+quote;
					*/
				}
			}
			
			if(dbin.number_tuples()<=0)
			{
				cerr<<"Waveform view has no data after joining: "<<
					"wfprocess+evlink+sclink+decon (if applicable)."<<endl;
				exit(-1);
			}
			
			if(use_arrival_data)
			{
				if(SEISPP_verbose) cout<<"Building catalog view ..."<<endl;
				DatascopeHandle ljhandle(dbin);
				try{
					ljhandle.lookup("event");
					ljhandle.natural_join("origin");
					ljhandle.natural_join("assoc");
					ljhandle.natural_join("arrival");
					//ljhandle.subset("sta==arrival.sta");
					ljhandle.subset("orid==prefor");
					ljhandle.subset("iphase=~/P/");
				}catch(SeisppError& serr)
				{
					cerr<<"Error in building: origin+assoc+arrival and subset with "
						<<"sta=arrival.sta & orid==prefor & phase=~/P/"<<endl;
					serr.log_error();
				}
				list<string> jk;
				jk.push_back("evid");
				jk.push_back("sta");
				dbin.join(ljhandle,jk,jk);
				if(SEISPP_verbose) cout<<"Number of rows after joining with catalog view: "
										<<dbin.number_tuples()<<endl;
				if(dbin.number_tuples()<=0)
				{
					cerr<<"Working view has no data after joining with: "
						<<"origin+assoc+arrival and subset with sta=arrival.sta "
						<<"& orid==prefor &phase=~/P/"<<endl;
					exit(-1);
				}
			}
        }
        			
		//read in netmag information.
		if(use_netmag_table)
		{
			try{
				if(SEISPP_verbose) cout<<"Joining with netmag table ..."<<endl;
				dbin.natural_join("netmag");
			}catch(SeisppError& serr)
			{
				cerr<<"Error in joining with netmag table!"<<endl;
				serr.log_error();
			}
		}
        
        //cerr<<"test"<<endl;
        if(dbin.number_tuples()<=0)
        {
            cerr << "No rows in input wfdisc/wfprocess table for database "
            	<<dbin_name<<endl;
            exit(-1);
        }
        logfile << "Number of rows in full wfdisc/wfprocess originally openned (may have duplicates) = "
            <<dbin.number_tuples()<<endl;
        /* Prep the input wfdisc table*/
        /* WARNING:  not adequate.  This needs to be changed as we 
           need evid.  May be able to fake this by searching for matching
           start times  */
        if(apply_subset)
        {
            cout << "Applying subset condition = [ "<<subset_condition<<" ]."<<endl;
            logfile << "Applied subset condition = [ "<<subset_condition<<" ]."<<endl;
            dbin.subset(subset_condition);
            cout << "Subset view number of rows = "<<dbin.number_tuples()<<endl;
        }
        list<string> sortkeys, groupkeys;
        sortkeys.push_back("sta");
        sortkeys.push_back("time");
        sortkeys.push_back("chan");
        groupkeys.push_back("sta");
        dbin.sort(sortkeys);
        dbin.group(groupkeys);
        cout << "Number of ensembles to process (grouped by sta)="
            << dbin.number_tuples()<<endl;
				/****************************************************************
				*****************************************************************
				*******<<<<<<<<<<<< PRE-FILTER SETUP >>>>>>>>>>********
				*****************************************************************
				*****************************************************************
				*/
		int nwsamp;
		double datadt, wavelet_width_parameter;
		WaveletNormalizationMethod nm=PEAK;
		TimeSeries wavelet;
		
        if(apply_prefilter) 
        {
        	logfile<<"Pre-filter type: "<<wavelet_type<<endl;
        	if(wavelet_type=="filter")
			{
				filterspec=control.get_string("filter");
				//debug
				if(SEISPP_verbose)
				{
					cout<<"Data will be pre-filered before stacking."<<endl;
					cout<<"Pre-filter specs = "<<filterspec<<endl;
				}
				logfile<<"Filter specs = "<<filterspec<<endl;
			}
			else if(wavelet_type=="gaussian")
			{
				nwsamp=control.get_int("wavelet_length");
				datadt=control.get_double("data_sample_interval");
				wavelet_width_parameter=control.get_double("wavelet_width_parameter");
				wavelet=gaussian_wavelet(nwsamp,datadt,
						wavelet_width_parameter,nm);
				//logfile<<"Pre-filter type: "<<wavelet_type<<endl;
			}
			// ricker wavelet is not suitable for RFeditor
			else if(wavelet_type=="ricker")
			{
				nwsamp=control.get_int("wavelet_length");
				datadt=control.get_double("data_sample_interval");
				wavelet_width_parameter=control.get_double("wavelet_width_parameter");
				wavelet=ricker_wavelet(nwsamp,datadt,
						wavelet_width_parameter,nm);
			}
			
			else
			{
				cout << "-"<<wavelet_type<<"-"<<endl;
				cout << "wavelet_type parameter="<<wavelet_type
					<< " is not supported."<<endl
					<< "Must be:  filter, gaussian"<<endl;
				usage();
			}
        	
        }
        		/****************************************************************
				*****************************************************************
				******<<<<<<<<<<<< BUILD GLOBAL EDITING OBJECTS >>>>>>>>>>*******
				*****************************************************************
				*****************************************************************
				*/
        /* This launches the editing windows */
        RFeditorEngine *rfe;
        //if(!GUIoff) 
        	rfe= new RFeditorEngine(trace_edit_params, GUIoff);
        /*TraceEditOperator object for trace editing.*/
        TraceEditOperator teo(trace_edit_params);
        int nsta=dbin.number_tuples(),ntrace(0),nradial(0);
        string sta,tracetype;
        		/****************************************************************
				*****************************************************************
				*************<<<<<<<<<<<< START MAIN LOOP >>>>>>>>>>*************
				*****************************************************************
				*****************************************************************
				*/
        vector<TimeSeries>::iterator im;
        TimeSeriesEnsemble radial,transverse,vertical,tse_edit0,tse_edit; //radial0
		FILE * fh_fa;
		if(get_FA)
		{
			fh_fa=fopen(FA_filename.c_str(),"w");
			if(use_arrival_data)
				fprintf(fh_fa,"STA    EVID    START_TIME    FA_LAG    FA_AMPR    SEAZ\n");
			else
				fprintf(fh_fa,"STA    EVID    START_TIME    FA_LAG    FA_AMPR\n");
		}
        for(i=0,dbin.rewind();i<nsta;++i,++dbin)
        {   
            cout <<">>++++++++++++++++++++++++++++++"<<endl
            	<<"Calling data reader for ensemble number ["<<i+1<<" / "<<nsta<<"]"<<endl;
            logfile<<">>++++++++++++++++++++++++++++++"<<endl
            	<<"Working on ensemble number ["<<i+1<<" / "<<nsta<<"]"<<endl;
            TimeSeriesEnsemble dall;
            ThreeComponentEnsemble dall_3c;//,dall_3c_bkp;
            /*
			*****************************************************************
			*****************************************************************
			<<<<<<<<<<<<< Preparing Working TimeSeries Ensemble >>>>>>>>>>>>>
			*****************************************************************
			*****************************************************************
			*/
            if(!datatype3c)  //TimeSeriesEnsemble input table
            {
            	TimeSeriesEnsemble dall0(dbin,mdl,mdlens,am);
            	dall=dall0;
            	dall0.member.clear();
            	sta=dall.get_string("sta");
            	ntrace=dall.member.size();
            	if(no_vertical_data) ntrace=ntrace/2;
            	else ntrace=ntrace/3;
            	tracetype=const_cast<char *>(tracetype_1c.c_str());
				logfile << "Read "<<dall.member.size()<<" "<<tracetype<<" RF traces for station = "
					<< sta <<endl;
				cout << "Read "<<dall.member.size()<<" "<<tracetype<<" RF traces for station = "
					<< sta <<endl;

				if(ntrace<minrfcutoff)
				{
					cout << "Station "<<sta<<" dropped.   Count below miminum of "
						<< minrfcutoff<<endl;
					logfile << "Station "<<sta
						<<" dropped. Count below miminum of "
						<< minrfcutoff<<endl;
					continue;
				}
				/* This routine sets a metadata item eventid based
				   on start time only.  Not a bulletproof approach but one 
				   that should work for RF data */
				if(SEISPP_verbose) cout<<"Setting event IDs ..."<<endl;
				set_eventids(dall);
				/* use start time as 0.  Also set moveout keyword
				to allow stacking in engine */
				//use user given data_start_time.
				double atime, t0, moveout;
				for(im=dall.member.begin();im!=dall.member.end();++im)
				{
					if(use_arrival_data)
        			{
						atime=im->get_double("atime");
						im->ator(atime-FA_reference_time);
						//DEBUG
						//cout<<"t0="<<im->t0<<", endtime="<<im->endtime()<<endl;
					}
					else
					{
						t0=im->t0;
						im->ator(t0);
					}
						//moveout=atime-t0;
						//cout<<"arrival time = "<<strtime(atime)<<endl;
						//cout<<im2->get_int("evid")<<"   "<<moveout<<endl; 
					
					im->put(moveout_keyword,0.0);
				}
				//TimeWindow twin0=teo.find_common_timewindow(dall);
                    //debug common window
                  //  cout<<"debug: twin.start="<<twin0.start<<", end="<<twin0.end<<endl;
                   // cout<<"debug: twin length="<<twin0.end-twin0.start<<endl;
				/* Could do this with the database, but I chose this 
				   algorithm because I think it will be more robust.
				   Main reason is I can use only a string fragment for
				   a match instead of demanding a full match */

				try 
				{
					radial=extract_by_chan(dall,rchan);
					cout << "Found "<<radial.member.size()
						<<" radial component RFs"<<endl;
					transverse=extract_by_chan(dall,tchan);
					cout << "Found "<<transverse.member.size()
						<<" transverse component RFs"<<endl;
					if(!no_vertical_data)
					{	if(save_vertical_channel || edit_on_vertical) 
						{
							vertical=extract_by_chan(dall,zchan);
							cout << "Found "<<vertical.member.size()
								<<" vertical component RFs"<<endl;
						}
					}
					else
					{
						vertical=radial;
						for(im=vertical.member.begin();im!=vertical.member.end();++im)
						{	
							im->put("chan",zchan);
							std::fill(im->s.begin(),im->s.end(),0.0);
						}
					}	
				}catch(SeisppError& serr)
				{
					cerr << "Problems in extract_by_chan.  Message "
						<<serr.what()<<endl
						<<"Skipping ensemble for station = "<<sta
						<<endl;
					continue;
				}
				
				dall.member.clear();	
				
				if(edit_on_radial) tse_edit=radial;
				else if(edit_on_transverse) tse_edit=transverse;
				else if(edit_on_vertical) tse_edit=vertical;
            }
            // load 3 component data from the view.
            //
            else
            {
            	//try{
            	ThreeComponentEnsemble dall0_3c(dbin, mdl,mdlens,am);
            	//}catch(SeisppError& serr)
            	//{
            	//	cerr<<"ERROR in builing ThreeComponentEnsemble:"<<endl;
            	//	serr.log_error();
            	//	exit(-1);
            	//}
            	dall_3c=dall0_3c;
            	dall0_3c.member.clear();
            	//if(!save_filtered_data && apply_prefilter) dall_3c_bkp=dall0_3c;
            	sta=dall_3c.get_string("sta");
            	ntrace=dall_3c.member.size();
            	//nradial=ntrace;
            	tracetype=const_cast<char *>(tracetype_3c.c_str());
            	logfile << "Read "<<ntrace<<" "<<tracetype<<" RF traces for station = "
					<< sta <<endl;
				cout << "Read "<<ntrace<<" "<<tracetype<<" RF traces for station = "
					<< sta <<endl;

				if(ntrace<minrfcutoff)
				{
					cout << "Station "<<sta<<" dropped.   Count below miminum of "
						<< minrfcutoff<<endl;
					logfile << "Station "<<sta
						<<" dropped. Count below miminum of "
						<< minrfcutoff<<endl;
					continue;
				}
				/* This routine sets a metadata item eventid based
				   on start time only.  Not a bulletproof approach but one 
				   that should work for RF data */
				if(SEISPP_verbose) cout<<"Setting event IDs ..."<<endl;
				set_eventids(dall_3c);
				/* use start time as 0.  Also set moveout keyword
				to allow stacking in engine */
				vector<ThreeComponentSeismogram>::iterator im2;
				double atime, t0, moveout;
				for(im2=dall_3c.member.begin();im2!=dall_3c.member.end();++im2)
				{
					//cout<<strtime(im2->t0)<<endl;
					if(use_arrival_data)
        			{
						atime=im2->get_double("atime");
						im2->ator(atime-FA_reference_time);
					}
					else
					{
						t0=im2->t0;
						im2->ator(t0);
					}
					//moveout=atime-t0;
					//cout<<"arrival time = "<<strtime(atime)<<endl;
					//cout<<im2->get_int("evid")<<"   "<<moveout<<endl; 
					
					im2->put(moveout_keyword,0.0);
				}
				
				//extract TimeSeries ensemble by chan codes.
				try
				{
					if(save_3C_data || review_mode)
					{
						shared_ptr<TimeSeriesEnsemble> data=ExtractComponent(dall_3c,edit_component);
						tse_edit=*data;
						data.reset();
						/*
						if(pre_edit_FA && !edit_on_radial) 
						{
							auto_ptr<TimeSeriesEnsemble> rdata=ExtractComponent(dall_3c,1);
							radial=*rdata;
							rdata.reset();
						}
						*/
					}
					else
					{
						shared_ptr<TimeSeriesEnsemble> tdata=ExtractComponent(dall_3c,0);
						transverse=*tdata;
						transverse.put("chan",tchan);
						shared_ptr<TimeSeriesEnsemble> rdata=ExtractComponent(dall_3c,1);
						radial=*rdata;
						radial.put("chan",rchan);
						tdata.reset();
						rdata.reset();
						if(save_vertical_channel) 
						{
							shared_ptr<TimeSeriesEnsemble> zdata=ExtractComponent(dall_3c,2);
							vertical=*zdata;
							vertical.put("chan",zchan);
							zdata.reset();
						}
						
						if(edit_on_radial) tse_edit=radial;
						else if(edit_on_transverse) tse_edit=transverse;
						else if(edit_on_vertical) tse_edit=vertical;
					}
				}catch(SeisppError& serr)
				{
					cerr << "Problems extracting TimeSeriesEnsemble by component.  Message "
						<<serr.what()<<endl
						<<"Skipping ensemble for station = "<<sta
						<<endl;
					continue;
				}				
            }
			
			tse_edit.put("chan",edit_chan_code);
			
            int j=set_duplicate_traces_to_false(tse_edit,false);
            if(j>0 && SEISPP_verbose)
            	cout<<"Duplicate traces in "<<edit_on_channel<<" (set to false) = "<<j<<endl;

            //kill timeseries with t0>0
            for(long i=0;i<tse_edit.member.size();i++)
            {
            	if(tse_edit.member[i].t0>FA_reference_time+MYZERO)
            		{
            			cout<<"***Set trace with t0 > FA_reference_time to FALSE!"<<endl
            				<<"    Start time:"
            				<<strtime(tse_edit.member[i].get_double(string("time")))<<endl;
            			tse_edit.member[i].live=false;
            			//transverse.member[i].live=false;
            			//if(pre_edit_FA && !edit_on_radial) radial.member[i].live=false;
            		}
            }
            set<long> evids_killed=teo.find_false_traces(tse_edit);
            //save original data before applying filters. Xiaotao Yang 1/22/2015
            if(SEISPP_verbose) cout<<"-- Excluding false traces ..."<<endl;
            /*
            if(pre_edit_FA && !edit_on_radial)
            {	radial0=teo.exclude_false_traces(radial);
            	radial=radial0;
            }*/
            //transverse0=teo.exclude_false_traces(transverse);
            //radial=radial0;
            //transverse=transverse0;
            tse_edit0=teo.exclude_false_traces(tse_edit);
            tse_edit=tse_edit0;
			cout<<"After excluding false traces: "<<edit_on_channel<<" = "<<tse_edit.member.size()<<endl;
			//	<<", transverse = "<<transverse.member.size()<<endl;
			// apply prefilter.Xiaotao Yang 01/12/2015
			cout<<"** Pre-filtering :: ";
			if(apply_prefilter) 
			{
				if(wavelet_type=="filter")
				{	
					TimeInvariantFilter filter(filterspec);
					logfile<<"Applying TimeInvariant Filter to "<<edit_on_channel<<" ensemble ..."<<endl;
					if(SEISPP_verbose) 
						cout<<"Applying TimeInvariant Filter to "<<edit_on_channel<<" ensemble ..."<<endl;
					//SEISPP::FilterEnsemble(radial,filter);
					SEISPP::FilterEnsemble(tse_edit,filter);
					//if(pre_edit_FA && !edit_on_radial) SEISPP::FilterEnsemble(radial0,filter);
				}
				else if(wavelet_type=="gaussian" || wavelet_type=="ricker")
				{
					//TimeWindow twin(radial.member[0].t0,radial.member[0].endtime());
                    TimeWindow twin=teo.find_common_timewindow(tse_edit);
                    //debug common window
                    cout<<"twin.start="<<twin.start<<", end="<<twin.end<<endl;
                    //teo.convolve_ensemble(wavelet,radial,true,&twin);
                    if(SEISPP_verbose) 
						cout<<"Convolving "<<edit_on_channel<<" ensemble with wavelet: "<<wavelet_type<<" ..."<<endl;
                    teo.convolve_ensemble(wavelet,tse_edit,true,&twin);
                    //if(pre_edit_FA && !edit_on_radial) teo.convolve_ensemble(wavelet,radial0,true,&twin);
				}
			}  // end of applying filter

			/*
			|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			<<<<<<<<<<<<<<<<<<<<<<<< Starting Editor >>>>>>>>>>>>>>>>>>>>>>>>
			|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
			*/
			if(get_FA && pre_edit_FA)
			//compute FA and write them out into a text file.
			{
				if(SEISPP_verbose) cout<<"Detecting first arrivals before applying editings ..."<<endl;
				TimeWindow FA_search_window=TimeWindow(FA_search_TW_start+FA_reference_time,
							FA_search_TW_end+FA_reference_time);
				vector<TimeSeries>::iterator iptr;
				TimeSeriesEnsemble tse_tmp(tse_edit);
				//if(!edit_on_radial) tse_tmp=radial0;
				for(iptr=tse_tmp.member.begin(); iptr!=tse_tmp.member.end();iptr++)
				{
					TimeSeries tmpts=teo.trim_data(*iptr,FA_search_window);
					string FA_type=find_FirstArrival(tmpts,FA_sensitivity,
								FA_detect_length,data_shaping_wavelet_type);
					double FA_time=tmpts.get_double(FA_time_key);
					double FA_lag=FA_time - FA_reference_time;
					double FA_amplitude=tmpts.get_double(FA_amplitude_key);
					int evid_tmp=iptr->get_int(evidkey);
					if(use_arrival_data)
						fprintf(fh_fa,"%5s   %6d    %15.3f   %6.3f   %8.4f   %6.1f\n",
							sta.c_str(),evid_tmp,iptr->get_double("time"),
							FA_lag,FA_amplitude,iptr->get_double(seaz_key));
					else
						fprintf(fh_fa,"%5s   %6d    %15.3f   %6.3f   %8.4f\n",
							sta.c_str(),evid_tmp,iptr->get_double("time"),
							FA_lag,FA_amplitude);
					
					tmpts.s.clear();
				}
				
				if(SEISPP_verbose) 
					cout<<"Saved [ "<<tse_tmp.member.size()
						<<" ] FA information for "<<edit_on_channel<<" data of station [ "<<sta<<" ]."<<endl;
				tse_tmp.member.clear();
				//exit(0);
			}
			
			if(!review_mode || !GUIoff)
			{
				if(SEISPP_verbose) cout<< "Loading data into editor ..."<<endl;
					//<< "When read, edit radial first then transverse."
					//<<endl;
		
				set<long> kills; //,tkills;
				if(!GUIoff)
					try{
						kills=rfe->edit(tse_edit); 
					}catch(...)
					{
						cerr<<"** Error in running GUI editor!"<<endl;
						exit(-1);
					}
				else
					try{
						kills=rfe->edit(tse_edit,trace_edit_params); 
					}catch(...)
					{
						cerr<<"** Error in running GUIoff editor!"<<endl;
						exit(-1);
					}
				rfe->save_statistics(logfilename);

				// stacktype for stack type. Xiaotao Yang
		
				/*
				=================================================================
				=================================================================
				<<<<<<<<<<<<<<<<<<<<<<<<<< Save or Not ? >>>>>>>>>>>>>>>>>>>>>>>>
				=================================================================
				=================================================================
				*/
				//find total number of killed traces.
				set<long> evids_killed2;
				evids_killed2=kills;
				//cout<<evids_killed2.size()<<endl;
				if(evids_killed2.size()>0)
					evids_killed.insert(evids_killed2.begin(),evids_killed2.end());
				cout<<"Found [ "<<evids_killed.size()<<" ] killed traces."<<endl;
				if(review_mode) 
				{
					cout<<"Review mode is on. Go to the next without saving the edits!"<<endl;
					rfe->reset_statistics();
				}
				else
				{
					//compute and save FA information
					if(get_FA && post_edit_FA )
					//compute FA and write them out into a text file.
					{
						if(evids_killed.size()>0)
						{
							if(SEISPP_verbose) cout<<"Getting FA: Applying kills to "
							<<edit_on_channel<<endl;
							teo.apply_kills(tse_edit,evids_killed);
						}
						TimeSeriesEnsemble tse_tmp=teo.exclude_false_traces(tse_edit);
						if(SEISPP_verbose) 
							cout<<"Getting FA: Detecting first arrivals after applying editings ..."<<endl;
						TimeWindow FA_search_window=TimeWindow(FA_search_TW_start+FA_reference_time,
									FA_search_TW_end+FA_reference_time);
						vector<TimeSeries>::iterator iptr;
						for(iptr=tse_tmp.member.begin(); iptr!=tse_tmp.member.end();iptr++)
						{
							TimeSeries tmpts=teo.trim_data(*iptr,FA_search_window);
							string FA_type=find_FirstArrival(tmpts,FA_sensitivity,
										FA_detect_length,data_shaping_wavelet_type);
							double FA_time=tmpts.get_double(FA_time_key);
							double FA_lag=FA_time - FA_reference_time;
							double FA_amplitude=tmpts.get_double(FA_amplitude_key);
							int evid_tmp=iptr->get_int(evidkey);
							if(use_arrival_data)
								fprintf(fh_fa,"%5s   %6d    %15.3f   %6.3f   %8.4f   %6.1f\n",
									sta.c_str(),evid_tmp,iptr->get_double("time"),
									FA_lag,FA_amplitude,iptr->get_double(seaz_key));
							else
								fprintf(fh_fa,"%5s   %6d    %15.3f   %6.3f   %8.4f\n",
									sta.c_str(),evid_tmp,iptr->get_double("time"),
									FA_lag,FA_amplitude);
							tmpts.s.clear();
						}
				
						if(SEISPP_verbose) 
							cout<<"Saved [ "<<tse_tmp.member.size()
								<<" ] FA information for "<<edit_on_channel
								<<" data of station [ "<<sta<<" ]."<<endl;
						tse_tmp.member.clear();
						//exit(0);
					}

					int nsaved;
					string outtable;
					////save edit statistics information first.
					if(save_edit_summary_to_file)
					{
						rfe->save_statistics_summary(editsummaryfname,csversion);
					}
					//save statistics to db table: tredit.
					rfe->save_statistics_summary(dbout,2,csversion);
					//reset is required to avoiding duplicate/accumulated statistical information.
					rfe->reset_statistics();
					//starting saving waveform data.
				
					if(save_3C_data)
					{
						if(SEISPP_verbose) 
							cout <<"Applying kills to 3C data ... "<<endl;
						teo.apply_kills(dall_3c, evids_killed);
						vector<ThreeComponentSeismogram>::iterator im2;
						if(apply_prefilter)
						{
							if(save_filtered_data)
							{
								if(wavelet_type=="filter")
								{
									TimeInvariantFilter filter(filterspec);
									SEISPP::FilterEnsemble(dall_3c,filter);
								}
								else if(wavelet_type=="gaussian" || wavelet_type=="ricker")
								{
									//TimeWindow twin(dall_3c.member[0].t0,dall_3c.member[0].endtime());
									TimeWindow twin=teo.find_common_timewindow(dall_3c);
									teo.convolve_ensemble(wavelet,dall_3c,true,&twin);
								}
							}
						}
						for(im2=dall_3c.member.begin();im2!=dall_3c.member.end();++im2)
						{	
							//change back to absolute time.
							if(use_arrival_data)
							{
								double atime=im2->get_double("atime");
								im2->rtoa(atime-FA_reference_time);
							}
							else
							{
								double t0=im2->get_double("time");
								im2->rtoa(t0);
							}
						}
						//save wfdisc table seperately if use wfprocess as input
                        if(save_wfdisc_table)
                        {
							if(SEISPP_verbose) cout <<"Saving to db (wfdisc). Please wait ..."<<endl;
							outtable="wfdisc";
							save_to_db(dall_3c,mdlout_wfd,
									am,dbout,outdir,outdfile_base,
									save_metadata_only,false,rchan,tchan,zchan,outtable);
						}
						//saving edits to db.
						if(SEISPP_verbose) cout <<"Saving to db (wfprocess). Please wait ..."<<endl;
								//Xiaotao Yang 1/16/2015
						outtable="wfprocess";
						nsaved=save_to_db(dall_3c,mdlout,
									am,dbout,outdir,outdfile_base,
									save_metadata_only,save_decon_table);
					}
					else
					{
						if(apply_prefilter)
						{
							if(!save_filtered_data)
							{
								//radial=radial0;
								//transverse=transverse0;
								// apply kills. moved from RFeditorEngine.cc to this place. Xiaotao Yang
								if(evids_killed.size()>0)
								{
									//debug
									if(SEISPP_verbose) cout<<"Applying kills to radial."<<endl;
									//cout<<"rkill size: "<<rkills.size()<<endl;
									teo.apply_kills(radial,evids_killed);
									
									if(SEISPP_verbose) cout<<"Applying kills to transverse."<<endl;
									teo.apply_kills(transverse,evids_killed);
									//kill vertical if turned on "save vertical channel"
									if(save_vertical_channel)
									{
										if(SEISPP_verbose) cout<<"Applying kills to vertical ..."<<endl;
										teo.apply_kills(vertical,evids_killed);
									}
								}
							}
							else  //save filtered data
							{
								if(wavelet_type=="filter")
								{	
									TimeInvariantFilter filter(filterspec);
									SEISPP::FilterEnsemble(radial,filter);
									SEISPP::FilterEnsemble(vertical,filter);
									SEISPP::FilterEnsemble(transverse,filter);
								}
								else if(wavelet_type=="gaussian" || wavelet_type=="ricker")
								{
									//TimeWindow twin(radial.member[0].t0,radial.member[0].endtime());
                                    TimeWindow twin=teo.find_common_timewindow(radial);
                                    teo.convolve_ensemble(wavelet,radial,true,&twin);
									teo.convolve_ensemble(wavelet,vertical,true,&twin);
									teo.convolve_ensemble(wavelet,transverse,true,&twin);
								}
								if(evids_killed.size()>0)
								{
									//debug
									//radial=radial0;
									if(SEISPP_verbose) cout<<"Applying kills to radial."<<endl;
									//cout<<"rkill size: "<<rkills.size()<<endl;
									teo.apply_kills(radial,evids_killed);
									
									if(SEISPP_verbose) cout<<"Applying kills to transverse."<<endl;
									teo.apply_kills(transverse,evids_killed);
									//kill vertical if turned on "save vertical channel"
									if(save_vertical_channel)
									{
										if(SEISPP_verbose) cout<<"Applying kills to vertical ..."<<endl;
										teo.apply_kills(vertical,evids_killed);
									}
								}
							}
						}
						else
						{	
							if(evids_killed.size()>0)
								{
									//debug
									if(SEISPP_verbose) cout<<"Applying kills to radial."<<endl;
									//cout<<"rkill size: "<<rkills.size()<<endl;
									teo.apply_kills(radial,evids_killed);
									
									if(SEISPP_verbose) cout<<"Applying kills to transverse."<<endl;
									teo.apply_kills(transverse,evids_killed);
									//kill vertical if turned on "save vertical channel"
									if(save_vertical_channel)
									{
										if(SEISPP_verbose) cout<<"Applying kills to vertical ..."<<endl;
										teo.apply_kills(vertical,evids_killed);
									}
								}
						}
						
						
						// change time reference to absolute.
						double t0;
						for(im=radial.member.begin();im!=radial.member.end();++im)
						{
							//im->rtoa(t0);
							if(use_arrival_data)
							{
								double atime=im->get_double("atime");
								im->rtoa(atime-FA_reference_time);
							}
							else
							{
								t0=im->get_double("time");
								im->rtoa(t0);
							}
							//Fragile way to handle this, but skip stack traces
							//if(!(im->live))
							//	logfile << "Deleting " <<sta<<":"<< im->get_string("chan")
							//		<< " for time "<<strtime(t0)<<endl;
						
						}
					
						for(im=transverse.member.begin();im!=transverse.member.end();++im)
						{
							//double t0=im->get_double("time");
							//im->rtoa(t0);
							if(use_arrival_data)
							{
								double atime=im->get_double("atime");
								im->rtoa(atime-FA_reference_time);
							}
							else
							{
								t0=im->get_double("time");
								im->rtoa(t0);
							}
						}
						if(save_vertical_channel)
							for(im=vertical.member.begin();im!=vertical.member.end();++im)
							{
								//double t0=im->get_double("time");
								//im->rtoa(t0);
								if(use_arrival_data)
								{
									double atime=im->get_double("atime");
									im->rtoa(atime-FA_reference_time);
								}
								else
								{
									t0=im->get_double("time");
									im->rtoa(t0);
								}
							}
						//saving edits to db.
						
						if(!use_wfdisc_in) 
						{	
							
							if(save_wfdisc_table)
							{
								if(SEISPP_verbose) cout <<"Saving to db (wfdisc). Please wait ..."<<endl;  
								nsaved=save_to_db(radial,transverse,vertical,mdlout_wfd,
										am,dbout,outdir,outdfile_base,
										save_metadata_only,save_vertical_channel,
										save_decon_table,rchan, tchan, zchan,"wfdisc");
							}
							outtable="wfprocess";
						}
						else
						{
							if(save_wfprocess_table)
							{
								if(SEISPP_verbose) cout <<"Saving to db (wfprocess). Please wait ..."<<endl;  
								nsaved=save_to_db(radial,transverse,vertical,mdlout_wfp,
										am,dbout,outdir,outdfile_base,
										save_metadata_only,save_vertical_channel,
										save_decon_table,rchan, tchan, zchan,"wfprocess");
							}
							outtable="wfdisc";
						}
						
						if(SEISPP_verbose) cout <<"Saving to db ("<<outtable<<"). Please wait ..."<<endl;
						nsaved=save_to_db(radial,transverse,vertical,mdlout,
										am,dbout,outdir,outdfile_base,
										save_metadata_only,save_vertical_channel,
										save_decon_table,rchan, tchan, zchan,outtable);
					}
					
					logfile << "Saved "<<nsaved<<" RFs for station "<<sta<<endl;
					if(SEISPP_verbose) cout << "Saved "<<nsaved<<" RFs for station "<<sta<<endl;
				}
				evids_killed2.clear();
				kills.clear(); 
            }
			//clear set containers.
			evids_killed.clear();   
			dall_3c.member.clear();    
        }
        radial.member.clear();
        transverse.member.clear();
        vertical.member.clear();
        tse_edit0.member.clear();
        tse_edit.member.clear();
        if(get_FA) fclose(fh_fa);
        if(SEISPP_verbose) cout<<"RFeditor finished."<<endl;
        logfile << "RFeditor finished on [ " << dbin_name << " ] "  <<endl;
        logfile.close();
        //delete rfe;
    }catch(SeisppError& serr)
    {
        serr.log_error();
    }
    catch(exception& stdexcept)
    {
        cerr << "Exception thrown:  "<<stdexcept.what()<<endl;
    }
}
//END OF RFEDITOR.
