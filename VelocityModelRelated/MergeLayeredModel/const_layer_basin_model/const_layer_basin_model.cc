#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "Metadata.h"
#include "gclgrid.h"
#include "LatLong-UTMconversion.h"
#include "GriddedBasinModel.h"
#include "seispp.h"
using namespace std;
GCLscalarfield *PetrelGrid2GCL(string fn,string utmzone,
        double lat, double lon, double r, double az);

GCLscalarfield *PetrelGrid2GCL(string fname,string utmzone,
        double lat0, double lon0, double r0, double az)
{
    const string base_error("PetrelGrid2GCL procedure:  ");
    ifstream inp;
    try{
        inp.open(fname.c_str(),ios::in);
        if(inp.fail())
        {
            throw SeisppError(
                string("PetrelExport2GCL:  Open failed on file")
                    + fname);
        }
        int n1,n2;
        char line[256];
        while(inp.getline(line,256))
        {
            string discard;
            if(line[0]!='#') break;
            string inputline(line);
            if(inputline.find("Grid_size:")!=std::string::npos)
            {
                stringstream ss(line);
                string dummy;
                ss >> dummy; //skip # character
                ss >> dummy; //skip Grid_size keyword
                ss >> n1;
                ss >> dummy; // skip x character
                ss >> n2;
            }
        }
        GCLscalarfield *result;
        try {
            /* This is kind of an ugly construct that is a relic.
               setting nominal spacings initially to 1.0 and 
               setting iorigin and jorigin to 0 - should be 
               harmless as they are relics. */
            GCLgrid grd(n1,n2,fname,lat0,lon0,r0,az,1.0,1.0,0,0);
            result=new GCLscalarfield(grd);
        }
        catch(GCLgridError &gerr)
        {
            // translate to SeisppError to avoid multiple handlers
            throw SeisppError(base_error
                    + "GCLgridError was thrown by constructors - message below\n"
                    + gerr.what());
        }

        int ellipsoid_number(23);  //For library routine used this is WGS-84
        double x,y,z,i,j;
        Geographic_point gp;
        int count(0);
        do {
            stringstream ss(line);
            ss >> x;
            ss >> y;
            ss >> z;
            ss >> i;  
            ss >> j;
            if( (i<1) || (i>n1) || (j<1) || (j>n2) )
            {
                throw SeisppError(base_error
                        + "input indices out of range - offending line below"
                        + line);
            }
            double latdeg,londeg;
            UTMtoLL(ellipsoid_number,y,x,utmzone.c_str(),latdeg,londeg);
            gp.lat=rad(latdeg);
            gp.lon=rad(londeg);
            gp.r=r0_ellipse(gp.lat);
            /* z values read are in m and with positive up.  Hence divide by
               1000 and add to r */
            gp.r += (z/1000.0);   
            Cartesian_point cp=result->gtoc(gp);
            int ii=i-1;
            int jj=j-1;  // indices input a fortran order - need C here
            result->x1[ii][jj]=cp.x1;
            result->x2[ii][jj]=cp.x2;
            result->x3[ii][jj]=cp.x3;
            //result->val[ii][jj]=z/1000.0;   // m conversion to km
            //cout<<"layer = "<<fname<<", x= "<<x<<",y= "<<y
            //	<<", lon="<<londeg<<", lat="<<latdeg<<",z/1000 = "<<result->val[ii][jj]<<endl;
            ++count;
        }while(inp.getline(line,256));
        if(count!=(n1*n2)) 
            throw SeisppError(base_error
                    + "File inconsistency - line count does not match header");
        return(result);
    }
    catch(...){throw;};
}

/* This small helper will eventually be put in seispp. Putting it here for
   a test drive.  Note may want to extend this with an optional
   arg for tag to use constructor with Arr tags for blocks of stuff */
Metadata BuildMetadataFromPf(string pffile) 
{
    Pf *pf;
    int iret;
    iret=pfread(const_cast<char *>(pffile.c_str()),&pf);
    if(iret!=0)
    {
        throw MetadataError(string("BuildMetadataFromPF procedure:  ")
                +"pfread failed returning nonzero error code");
    }
    Metadata result(pf);
    return(result);
}
void usage()
{
    cerr << "const_layer_model modelname [-pf pffile] < listfile"<<endl
        << " infile item order:  isopach_file vp vs rho"<<endl;
    exit(-1);
}
bool SEISPP::SEISPP_verbose(true);
int main(int argc, char **argv)
{
    string pfname("const_layer_basin_model");
    if(argc<2) usage();
    /* used for file name and to give a tag to model */
    string modelname(argv[1]);
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
        else
            usage();
    }
    try {
        Metadata control;
        control=BuildMetadataFromPf(pfname);
        double lat0,lon0,r0,az;
        lat0=control.get_double("origin_latitude");
        lon0=control.get_double("origin_longitude");
        lat0=rad(lat0);
        lon0=rad(lon0);
        r0=control.get_double("origin_radius");
        az=control.get_double("azimuth_y_axis");
        az=rad(az);
        string utmzone=control.get_string("utmzone");
        string outdir=control.get_string("output_directory");
        vector<GCLscalarfield> isopachs;
        vector<double> vp,vs,rho;
        string isofile;
        double vp_layer,vs_layer,rho_layer;
        char line[256];
        while(cin.getline(line,256))
        {
            GCLscalarfield *current;
            stringstream ss(line);
            ss >> isofile;
            ss >> vp_layer;
            ss >> vs_layer;
            ss >> rho_layer;
            current = PetrelGrid2GCL(isofile,utmzone,lat0,lon0,r0,az);
            isopachs.push_back(*current);
            vp.push_back(vp_layer);
            vs.push_back(vs_layer);
            rho.push_back(rho_layer);
            delete current;
        }
        GriddedBasinModel model(isopachs,vp,vs,rho,modelname);
        /* Save method is just the GCLvectorfield3d method */
        model.save(modelname,outdir);
        //DEBUG
        
        cout << "Geometry check "<<endl;
        double zsum;
        for(int i=0;i<model.n1;++i)
            for(int j=0;j<model.n2;++j)
            {
                zsum=0.0;
                for(int k=model.n3-1;k>0;--k)
                {
                    zsum += (model.depth(i,j,k-1)) - (model.depth(i,j,k));
                    cout <<model.depth(i,j,k)<<" "<<zsum<<" ";
                }
                cout <<endl;
                cout <<"comparison  "<< i<<" "<<j<<" "<<zsum<<" "<<model.depth(i,j,0)
                    <<" "<<model.depth(i,j,0)-zsum<<endl;
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
}
