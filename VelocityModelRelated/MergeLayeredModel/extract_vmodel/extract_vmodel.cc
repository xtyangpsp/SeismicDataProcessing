#include <string>
#include <fstream>
#include "perf.h"
#include "gclgrid.h"
#include "seispp.h"
using namespace std;
using namespace SEISPP;

/*
Modified from original ILBasinCrust.cc by Gary Pavlis.
The model file must be a GCLvectorfield3d grid with four components for vp, vs, density and thickness.
This program extracts 1-d layered model for each point (usually station) with output model formated as [migsimulation] input for StaVariableLayeredSynthetic method.

Xiaotao Yang @IU
2016.2.09

2016.5.11 XT Yang @ IU
	Added option of outputing to other files instead of stdout. In this case, the program generate a log file containing status of lookup for each station (0: failed; 1: found).

*/
void usage()
{
    cerr << "extract_vmodel gclgridmodelfile [-tp|--top-pad pad_length -o outfile -v] < pointlistfile"<<endl
        << " * pointlistfile contains: sta  lat lon elev(km) - normally from site table"<<endl
        << " * if use other output file, output the status to sdout: sta 0(failed)|1(found)"<<endl
        << " * -tp|--top-pad: apply top pad with specified pad_length, default 5.0."<<endl;
    exit(-1);
}
void pad_top_model(GCLvectorfield3d& f, double pad_length)
{
    //const double pad_length(5.0);
    Cartesian_point cp;
    Geographic_point gp;
    int k0=f.n3 - 1;
    int i,j;
    for(i=0;i<f.n1;++i)
        for(j=0;j<f.n2;++j)
        {
            gp=f.geo_coordinates(i,j,k0);
            gp.r+=pad_length;
            cp=f.gtoc(gp);
            f.x1[i][j][k0]=cp.x1;
            f.x2[i][j][k0]=cp.x2;
            f.x3[i][j][k0]=cp.x3;
        }
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    if(argc<2) usage();
    string modelgridfile(argv[1]);
    string outfile("out.txt");
    double pad_length(5.0);
    bool apply_top_pad(false);
    ofstream outstrm;
	bool out_to_other(false);
    int i;
    for(i=2;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-tp" || sarg=="--top-pad")
        {
            ++i;
            if(i>=argc) pad_length=5.0;
            else pad_length=atof(argv[i]);
            apply_top_pad=true;
        }
        else if(sarg=="-o")
		{
			++i;
			if(i>=argc) usage();
			out_to_other=true;
			outfile=string(argv[i]);
		}
		else if(sarg=="-v")
			SEISPP_verbose=true;
        else
            usage();
    }
    
    if(out_to_other)
    {
    	try
    	{outstrm.open(outfile.c_str(),ios::out);}catch (ios::failure& var)
			{
				cerr << "Open failure on outfile="<<outfile
					<<endl;
					
				if(SEISPP_verbose) cerr<< "System message:  "
					<< var.what() <<endl;
				usage();
			}
    }
    try {
        int k,kk;
        //Crust1_0 c1p0;
        
        GCLvectorfield3d vmodel(modelgridfile);
        /* This procedure is super bubba, but it is necessary because 
           the upper layers in this model are very very thin and the 
           lookup method used in the gclgrid library is not appropriate
           for this geometry.   */
        if(apply_top_pad) pad_top_model(vmodel,pad_length);
        /*storing the complete model in this matrix.  We 
        append 4 lower layers from crust1.0 to basin model.
        4 columns for vp,vs,rho, and thickness. -1 is explicitly
        used because the last surface in vmodel is basement 
        that replaces comparable surface in crust1.0 */
        int outsize=vmodel.n3;
        dmatrix modeldata(outsize,4);
        char inpline[128];
        while(cin.getline(inpline,128))
        {
            stringstream ss(inpline);
            string sta;
            double lat,lon,latdeg,londeg; //,elev;
            ss >> sta;
            ss >> latdeg;
            ss >> londeg;
            //ss >> elev;
            /* Find this position in the GCLgrid and extract the 
               vector of data for the basin model for each layer.
               A bit odd because we only use the lookup function
               for the first point */
            lat=rad(latdeg);
            lon=rad(londeg);
            //elev /= 1000.0;
            double r=r0_ellipse(lat);
            Cartesian_point cp=vmodel.gtoc(lat,lon,r);
            //bool use_only_crust1p0(false);
            if(vmodel.lookup(cp.x1,cp.x2,cp.x3))
            {
                if(SEISPP_verbose) cerr << "Lookup failed for sta "<<sta<<endl
                    << "Skipped this station"<<endl;
                //use_only_crust1p0=true;
                if(out_to_other) cout<<sta<<"    "<<0<<endl;
                continue;
            }
            else
            	if(out_to_other) cout<<sta<<"    "<<1<<endl;
            
            int position[3];
            //if(!use_only_crust1p0)
            //{
            vmodel.get_index(position);
            /* now copy the vmodel values to the output matrix*/
            int kk;
            for(kk=0,k=vmodel.n3-1;k>=0;++kk,--k)
            {
                for(int l=0;l<4;++l)
                    modeldata(kk,l)=vmodel.val[position[0]][position[1]][k][l];
            }
            if(out_to_other) 
            {
            	outstrm << sta << " " << outsize<<endl;
            	outstrm << modeldata;
            }
            else
            {
            	cout << sta << " " << outsize<<endl;
            	cout << modeldata;
			}
            //DEBUG
        }
    }catch(SeisppError& serr)
    {
        serr.log_error();
        exit(-1);
    }
    catch(std::exception& err)
    {
        cerr << err.what()<<endl;
        exit(-1);
    }
    catch(...)
    {
        cerr << "Something threw an exception I do not know how to handle"<<endl;
        exit(-1);
    }
}
