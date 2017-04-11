#include <vector>
#include <memory>
#include "stock.h"
#include "pf.h"
#include "seispp.h"
#include "Metadata.h"
#include "EventCatalog.h"
#include "dbpp.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
void usage()
{
	cerr << "telescluster db [-pf pffile -V]" <<endl;
	exit(-1);
}
/* Simple radial grid object used interally here */
class RadialGrid
{
public:
	double lat0,lon0;
	int naz,ndelta;
	int nazbins, ndelbins;
	vector<double> delta;
	vector<double> azimuth;

	RadialGrid(Pf *pf);
	/* note constructor lat0,lon0 should be passed as 
	degrees. Converted to radians */
	RadialGrid(double azmin, double azmax, int nazin,
		double delmin, double delmax, int ndelin,
		double lat0in, double lon0in);
	/* Recommended way to get a point.  Returns geographic coordinate location of a point */
	Geographic_point grid_point(int ir, int id);
	/* Return latitude in radians of a point */
	double lat(int ir, int id);
	/* same for longitude */
	double lon(int ir, int id);
};
RadialGrid::RadialGrid(Pf *pf)
{
	string base_error("RadialGrid constructor:  ");
	lat0=pfget_double(pf,"origin_latitude");	
	lon0=pfget_double(pf,"origin_longitude");	
	lat0=rad(lat0);
	lon0=rad(lon0);
	Tbl *t;
	t=pfget_tbl(pf,"delta_grid_points");
	if(t==NULL) throw SeisppError(base_error
		+string("delta_grid_points parameter is not in parameter file"));
	char *s;
	double valin;
	int i;
	for(i=0;i<maxtbl(t);++i)
	{
		s=(char *)gettbl(t,i);
		stringstream ss(s);
		ss >> valin;
		if( (valin<0.0) || (valin>180.0) )
			throw SeisppError(base_error
				+ string("Illegal distance specification=")
				+ s);
		delta.push_back(rad(valin));
	}
	ndelta=delta.size();
	ndelbins=ndelta-1;
	freetbl(t,0);
	t=pfget_tbl(pf,"azimuth_grid_points");
	if(t==NULL) throw SeisppError(base_error
			+ string("azimuth_grid_points parameter is not in parameter file"));
	for(i=0;i<maxtbl(t);++i)
	{
		s=(char *)gettbl(t,i);
		istringstream ss(s);
		ss >> valin;
		if( (valin<-180.0) || (valin>360.0) )
			throw SeisppError(base_error
			 + string("Illegal azimuth specification =")
			 + s);
		azimuth.push_back(valin);
	}
	naz=azimuth.size();
	nazbins=naz-1;
}
RadialGrid::RadialGrid(double azmin, double azmax, int nazin,
		double delmin, double delmax, int ndelin,
		double lat0in, double lon0in)
{
	const string base_error("RadialGrid parameterized constructor:  ");
	if(azmin>=azmax) throw SeisppError(base_error
		+ "Illegal azimuth range specified");
	if(delmin>=delmax) throw SeisppError(base_error
		+ "Illegal distance range specified");
	if( (nazin<=0) || (ndelin<=0) ) throw SeisppError(base_error
		+ "Number of points in azimuth and delta must be positive");
	naz=nazin;
	ndelta=ndelin;
	nazbins=naz-1;
	ndelbins=ndelta-1;
	lat0=rad(lat0in);
	lon0=rad(lon0in);
	double daz=(azmax-azmin)/static_cast<double>(nazbins);
	double ddel=(delmax-delmin)/static_cast<double>(ndelbins);
	int i;
	for(i=0;i<naz;++i) 
	  azimuth.push_back(rad(azmin+daz*static_cast<double>(i)));
	for(i=0;i<ndelta;++i) 
	  delta.push_back(rad(delmin+ddel*static_cast<double>(i)));
}
double RadialGrid::lat(int ir, int id)
{
    try {
	Geographic_point p=this->grid_point(ir,id);
	return(p.lat);
    } catch (...) {throw;};
}
double RadialGrid::lon(int ir, int id)
{
    try {
	Geographic_point p=this->grid_point(ir,id);
	return(p.lon);
    } catch (...) {throw;};
}
Geographic_point RadialGrid::grid_point(int ir, int id)
{
	if( (ir<0) || (ir>=nazbins) || (id<0) || (id>=ndelbins) )
	{
		stringstream ss;
		ss << "RadialGrid:  requested index (az,del)=("
			<<ir<<", "<<id<<") out of range"<<endl
			<<"Allowed: distance=(0,"<<ndelbins
			<<") azimuth=(0,"<<nazbins<<")"<<endl;
		throw SeisppError(ss.str());
	}
	double plat,plon,pr;
	double delcenter,azcenter;
	delcenter=(delta[id]+delta[id+1])/2.0;
	azcenter=(azimuth[ir]+azimuth[ir+1])/2.0;
	latlon(lat0,lon0,delta[id],azimuth[ir],&plat,&plon);
	latlon(lat0,lon0,delcenter,azcenter,&plat,&plon);
	pr=r0_ellipse(plon);
	Geographic_point result;
	result.lat=plat;
	result.lon=plon;
	result.r=pr;
	return(result);
}
/* This function object is used as the predicate for the subset method template in 
EventCatalog.  */
class SectorTest : public unary_function<Hypocenter,bool>
{
public:
	bool operator() (const Hypocenter& h) const
	{
		double del,az;
		dist(lat0,lon0,h.lat,h.lon,&del,&az);
		if( (del<delmin) || (del>delmax) ) return false;
		if( (az<azmin) || (az>azmax) ) return false;
		return true;
	};
	SectorTest(RadialGrid& g,int ia, int id);
private:
	double delmin,delmax;
	double azmin,azmax;
	double lat0,lon0;
};
SectorTest::SectorTest(RadialGrid& g, int ia, int id)
{
	const string base_error("SectorTest constructor:  ");
	if( (ia<0) || (ia>=g.nazbins) ) 
		throw SeisppError(base_error
		 + "Requested azimuth index outside allowed range");
	if( (id<0) || (id>=g.ndelbins) ) 
		throw SeisppError(base_error
		 + "Requested distance index outside allowed range");
	lat0=g.lat0;
	lon0=g.lon0;
	azmin=g.azimuth[ia];
	delmin=g.delta[id];
	azmax=g.azimuth[ia+1];
	delmax=g.delta[id+1];
}

/* Saves cluster table.  Returns gridid as nextid value that defines gridid */
long save_cluster_data(EventCatalog& events, DatascopeHandle& dbh, string gridname,long grididin)
{
	long gridid,evid;
	Metadata aux;
        /* This is used to signal to use dbnextid.  A potentially confusing hack fix made to 
           allow external setting of gridid */
        if(grididin<0)
	    gridid=dbnextid(dbh.db,"gridid");
        else
            gridid=grididin;
	int i;
	events.rewind();
	for(i=0;i<events.size();++i)
	{
		aux=events.current_aux();
		evid=aux.get_int("evid");
		dbh.append();
		dbh.put("gridname",gridname);
		dbh.put("gridid",gridid);
		dbh.put("evid",evid);
		++events;
	}
	return(gridid);
}
/* This routine both computes and saves a hypocentroid for a radial grid cell.  
lat0 and lon0 are assumed to be a reference cell grid point on one of the edges
of the cell.  We compute the hypocentroid using effectively easting and northing
distances from lat0 and lon0. This will be inaccurate if lat0 and lon0 are far
from the contents of the EventCatalog events object*/
void save_hypocentroid_data(EventCatalog& events, DatascopeHandle& dbh,
	double lat0, double lon0,string gridname, int gridid)
{
	int nevents=events.size();
	double deltaN(0.0), deltaE(0.0);
	double zmean(0.0);
	double del,az;
	events.rewind();
	int i;
	for(i=0;i<nevents;++i)
	{
		Hypocenter h=events.current();
		dist(lat0,lon0,h.lat,h.lon,&del,&az);
		deltaN+=del*cos(az);
		deltaE+=del*sin(az);
		zmean+=h.z;
		++events;
	}
	deltaN /= static_cast<double>(nevents);
	deltaE /= static_cast<double>(nevents);
	zmean /= static_cast<double>(nevents);
	del=sqrt(deltaN*deltaN + deltaE*deltaE);
	az=atan2(deltaE,deltaN);
	double hclat,hclon;
	latlon(lat0,lon0,del,az,&hclat,&hclon);
	dbh.append();
	dbh.put("gridname",gridname);
	dbh.put("gridid",gridid);
	dbh.put("hclat",deg(hclat));
	dbh.put("hclon",deg(hclon));
	dbh.put("hcdepth",zmean);
	dbh.put("dlat",deg(lat0));
	dbh.put("dlon",deg(lon0));
	dbh.put("depth",0.0);
	dbh.put("nass",nevents);
}

bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	string pfname("telecluster");
	if(argc<2) usage();
	string dbname(argv[1]);
	int i;
	for(i=2;i<argc;++i)
	{
		string sarg(argv[i]);
		if(sarg=="-pf")
		{
			++i;
			if(i>=argc) usage();
			pfname=string(argv[i]);
		}
		else if(sarg=="-V")
			SEISPP_verbose=true;
		else
			usage();
	}
	if(SEISPP_verbose) cout << "Program "<<argv[0] <<" begins processing db="<<dbname<<endl;
	Pf *pf;
	if(pfread(const_cast<char *>(pfname.c_str()),&pf))
	{
		cerr << "pfread failed on file="<<pfname<<endl;
		usage();
	}
	try {
		Metadata control(pf);
		/* Parse control parameters first */
		double origin_lat=control.get_double("origin_latitude");
		double origin_lon=control.get_double("origin_longitude");
		RadialGrid *grid;
		bool use_regular_grid=control.get_bool("use_regular_grid");
		string gridname=control.get_string("gridname");
		if(use_regular_grid)
		{
			double azmin=control.get_double("grid_minimum_azimuth");
			double azmax=control.get_double("grid_maximum_azimuth");
			int naz=control.get_int("number_grid_points_for_azimuth");
			double delmin=control.get_double("grid_minimum_delta");
			double delmax=control.get_double("grid_maximum_delta");
			int ndel=control.get_int("number_grid_points_for_delta");
			grid=new RadialGrid(azmin,azmax,naz,delmin,delmax,ndel,
				origin_lat,origin_lon);
		}
		else
		{
			grid=new RadialGrid(pf);
		}
		origin_lat=rad(origin_lat);  origin_lon=rad(origin_lon);
		string ttmethod=control.get_string("ttmethod");
		string ttmodel=control.get_string("ttmodel");
                /* Added option Nov 2012 for geographic merge capability.
                   When true gridid values will be computed from grid geometry
                   instead of using dbnextid. */
                bool use_count_for_gridid=control.get_bool("use_count_for_gridid");
		DatascopeHandle dbh(dbname,false);
		DatascopeHandle dbhcluster(dbh);
		dbhcluster.lookup("cluster");
		DatascopeHandle dbhypoc(dbh);
		dbhypoc.lookup("hypocentroid");
		/* This view is needed for EventCatalog constructor */
		dbh.lookup("event");
		dbh.natural_join("origin");
		dbh.subset("orid==prefor");
		/* Hand code this one required aux attribute */
		Metadata_typedef mdevid={string("evid"),MDint};
		MetadataList mdl;
		mdl.push_back(mdevid);
		AttributeMap am("css3.0");
		EventCatalog catalog(dynamic_cast<DatabaseHandle&>(dbh),
			mdl,am,ttmethod,ttmodel);
		if(SEISPP_verbose) cout << "Loaded event catalog with "
				<< catalog.size()
				<< " events"<<endl;
		long gridid;
		double latp,lonp;
		/* Main loop  over azimuth and delta segments */
		for(int ia=0;ia<(grid->naz)-1;++ia)
		    for(int id=0;id<(grid->ndelta)-1;++id)
		    {
		    	SectorTest thissector(*grid,ia,id);
		    	auto_ptr<EventCatalog> evsubset=catalog.subset<SectorTest>(thissector);
			if(evsubset->size() > 0)
			{
				Geographic_point p=grid->grid_point(ia,id);
				if(SEISPP_verbose) 
				  cout <<"Grid point at lat="
					<< deg(p.lat)
					<< " lon="
					<< deg(p.lon) 
					<< " has " << evsubset->size()
					<< " events"<<endl;
                                long test_gridid;
                                /* This is not an elegant way to handle this, but was
                                   done as a patch */
                                if(use_count_for_gridid)
                                    gridid=ia*(grid->naz)+id+1;
                                else
                                    gridid=-1;
				test_gridid=save_cluster_data(*evsubset,dbhcluster,gridname,gridid);
                                if(!use_count_for_gridid)  gridid=test_gridid;
				save_hypocentroid_data(*evsubset,dbhypoc,
					p.lat,p.lon,gridname,gridid);
			}
		    }
	} catch (SeisppError& serr)
	{
		serr.log_error();
		usage();
	}
        catch (...)
        {
            cerr << "Unknown exception thrown"<<endl;
        }
}
