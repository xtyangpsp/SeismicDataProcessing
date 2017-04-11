#include <vector>
#include "seispp.h"
#include "dbpp.h"
#include "gclgrid.h"
#include "VelocityModel_1d.h"
using namespace std;
using namespace SEISPP;

void usage()
{
	cerr << "project1dmod db gridname vmodfile "
	<< "[-field fieldname -mt modtype -p rayparameter -fileoutput -V]"<<endl
	<< "Default fieldname=vmodfile, modtype=P, rayparameter=0(s/km)"<<endl
        << "-fileoutput writes to a file with base name fieldname instead of db"<<endl;;
	exit(-1);
}
bool SEISPP::SEISPP_verbose(false);
int main(int argc, char **argv)
{
	if(argc<4) usage();
	string dbname(argv[1]);
	string gridname(argv[2]);
	string vmodname(argv[3]);
	string modtype("P");
	string fieldname=vmodname;
        bool fileoutput(false);
	double rayp(0.0);
	int i;
	for(i=4;i<argc;++i)
	{
		string argstr(argv[i]);
		if(argstr=="-field")
		{
			++i;
			fieldname=string(argv[i]);
		}
		else if(argstr=="-mt")
		{
			++i;
			modtype=string(argv[i]);
		}
		else if(argstr=="-p")
		{
			++i;
			rayp=atof(argv[i]);
		}
                else if(argstr=="-fileoutput")
                {
                        fileoutput=true;
                }
		else if(argstr=="-V")
		{
			SEISPP_verbose=true;
		}
		else
		{
			cerr << "Illegal argument="<<argstr<<endl;
			usage();
		}
	}
	try {
		cout << "Building database handles"<<endl;
		DatascopeHandle dbh(dbname,false);
		DatascopeHandle dbhgrid(dbh);
		dbhgrid.lookup("gclgdisk");
		cout << "Loading gridname="<<gridname<<endl;
		GCLgrid3d grid(dbhgrid,gridname);
		GCLscalarfield3d mod(grid);
		cout << "Loading 1d velocity model with file name="<<vmodname<<endl;
		VelocityModel_1d Vmod(vmodname,string("mod1d"),modtype);
		/* the following works only when the grid has constant
		depth slices for index 3 constant */
		vector<double> zi,vi;
		cout << "1D velocities (z,v,vapparent)"<<endl;
		for(i=mod.n3-1;i>=0;--i)
		{
			double z;
			z=mod.depth(0,0,i);
			zi.push_back(z);
			double vel,slow;
			vel=Vmod.getv(z);
			cout <<z<<" "<< vel <<"  ";
			if(rayp>0.0)
			{
				slow=1/vel;
				if(slow<rayp)
				{
					cerr << "Ray parameter mismatch"<<endl
					<< "slowness="<<slow<<" at z="<<z
					<< " but rayp="<<rayp<<endl
					<< "ray parameter must be less than slowness"<<endl;
					exit(-1);
				}
				slow=sqrt(slow*slow-rayp*rayp);
				vel=1.0/slow;
			}
			cout <<vel<<endl;
			vi.push_back(vel);
		}
		initialize_1Dscalar(mod,vi,zi);
		if(SEISPP_verbose)
		{
			cout << "Model written to db"<<endl;
			cout << mod;
		}
		/* for now this is a frozen directory name */
		string fielddir("vmodels");
                if(fileoutput)
                {
                    /* output name is set to field name */
                    mod.save(fieldname,fielddir);
                }
                else
                {
                    // This is db output mode.   
		    string gclgdir("");  /* null as save is not needed*/
                    /* We use the fieldname as the file name.   A bit
                       restrictive, but better than adding the baggage of a
                       pf file */
		    mod.save(dbhgrid,gclgdir,fielddir,fieldname,fieldname);
                }
	} catch (SeisppError& serr)
	{
		serr.log_error();
	}
	catch (std::exception& sterr)
	{
		cerr << sterr.what();
	}
}
