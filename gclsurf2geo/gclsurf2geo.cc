#include <fstream>
#include <iostream>
#include <sstream>
#include "stock.h"
#include "dbpp.h"
#include "TimeSeries.h"
#include "Metadata.h"
#include "seispp.h"
#include "gclgrid.h"
using namespace std;
using namespace SEISPP;
/* This is modified from surf_vtk_converter (by Yinzhi Wang), used to convert exported kingdom surface file into points in vtk format.
	The format of input should be inline crossline depth
   */
   
/* This modified version is designed to convert surface data in gclgrid format (j,i,attribute) to geographic points with the format of (lon, lat, attribute).

Changes:
1. program name is changed from surf_vtk_converter to gclsurf2geo
2. gclgrid information is read in from grid file directly without need of Antelope database.
*/

void usage()
{
	cerr << "gclsurf2geo gridfile [-v] < infile"<<endl
	      << "* output is stdout"<<endl
	      << "* gridfile shoule be 3d grid basename without extension."<<endl;

	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	int i,j;
	ios::sync_with_stdio();
	if(argc<2) usage();
	//string dbname(argv[1]);
	string gridfile(argv[1]);
	//string fieldname(argv[3]);
	ofstream outstrm;
	for(i=2;i<argc;++i)
	{
		string argstr=string(argv[i]);
		if(argstr=="-v")
			SEISPP_verbose=true;
		else
		{
			cerr << "Unknown argument = "<<argstr<<endl;
			usage();
		}
	}
	try 
	{
        /*
        DatascopeHandle dbh(dbname,true);
        dbh.lookup(string("gclgdisk"));
	    DatascopeHandle dbhg(dbh);
	    dbhg.lookup(string("gclfield"));
	    Dbptr db=dbh.db;
	    Dbptr dbgrd=dbhg.db;
	   
	    GCLvectorfield3d g(dbh,gridname,fieldname,5);
	     */
	    GCLgrid3d g(gridfile);
	    //GCLvectorfield3d g(gridfile);
        /* Now loop over grid.  */
        int i,j;
        double xi,xj,depth,lat,lon;
        double x1max,x2max,dxi,dxj;
        string temp_line,line;
        const double align_tolerance(0.01);
        bool grid_aligned(false);
        x1max=static_cast<double>(g.n1)-1;
        x2max=static_cast<double>(g.n2)-1;
        int nl(1);
        while(getline(cin,temp_line))
        {
         	if(SEISPP_verbose) cerr<<"Working on line number: "<<nl<<endl;
         	stringstream line(temp_line);
           	line>>xj;
           	line>>xi;
           	line>>depth;
    		i=static_cast<int>(xi);
    		j=static_cast<int>(xj);
    		dxi=fabs(xi-static_cast<double>(i));
    		dxj=fabs(xj-static_cast<double>(j));
    		if(dxi<align_tolerance && dxj<align_tolerance)
    			grid_aligned=true;
    		if(grid_aligned)
    		{
    			if(SEISPP_verbose) cerr<<"Grid is aligned. "<<endl;
    			lat=g.lat(i,j,g.n3-1);
    			lon=g.lon(i,j,g.n3-1);
    			lat=deg(lat);   
    			lon=deg(lon);
    		}
    		else
    		{
	    		if(SEISPP_verbose) cerr<<"Grid is NOT aligned. "<<endl;
	    		Cartesian_point cp;
    			Geographic_point gp=g.geo_coordinates (i,j,g.n3-1);
    			Cartesian_point cp1=g.gtoc(gp);
    			gp=g.geo_coordinates (i+1,j,g.n3-1);
    			Cartesian_point cp2=g.gtoc(gp);
    			gp=g.geo_coordinates (i,j+1,g.n3-1);
    			Cartesian_point cp3=g.gtoc(gp);
    			cp.x1=(cp2.x1-cp1.x1)*dxi+(cp3.x1-cp1.x1)*dxj+cp1.x1;
    			cp.x2=(cp2.x2-cp1.x2)*dxi+(cp3.x2-cp1.x2)*dxj+cp1.x2;
    			cp.x3=(cp2.x3-cp1.x3)*dxi+(cp3.x3-cp1.x3)*dxj+cp1.x3;
    			gp=g.ctog(cp);
    			lat=gp.lat;
    			lon=gp.lon;
    			lat=deg(lat);   
    			lon=deg(lon);
    		}
    		cout<<lon<<"	"<<lat<<"	"<<depth<<endl;
    		++nl;
        }
	}
	catch (SeisppError& serr)
	{
		serr.log_error();
	}
}
