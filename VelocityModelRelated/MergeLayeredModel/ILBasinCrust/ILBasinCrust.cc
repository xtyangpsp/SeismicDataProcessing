#include <string>
#include <fstream>
#include "perf.h"
#include "gclgrid.h"
#include "Crust1_0.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "ILBasinCrust GCLbasinfile < infile"<<endl
        << "  infile contains: sta  lat lon elev(m) - normally from site table"<<endl;
    exit(-1);
}
void pad_top_model(GCLvectorfield3d& f)
{
    const double pad_length(5.0);
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
    if(argc!=2) usage();
    try {
        int k,kk;
        Crust1_0 c1p0;
        string basinfile(argv[1]);
        GCLvectorfield3d basinmodel(basinfile);
        /* This procedure is super bubba, but it is necessary because 
           the upper layers in this model are very very thin and the 
           lookup method used in the gclgrid library is not appropriate
           for this geometry.   */
        pad_top_model(basinmodel);
        /*storing the complete model in this matrix.  We 
        append 4 lower layers from crust1.0 to basin model.
        4 columns for vp,vs,rho, and thickness. -1 is explicitly
        used because the last surface in basinmodel is basement 
        that replaces comparable surface in crust1.0 */
        int outsize=basinmodel.n3+4-1;
        dmatrix mergedmodel(outsize,4);
        char inpline[128];
        while(cin.getline(inpline,128))
        {
            stringstream ss(inpline);
            string sta;
            double lat,lon,latdeg,londeg,elev;
            ss >> sta;
            ss >> latdeg;
            ss >> londeg;
            ss >> elev;
            /* Find this position in the GCLgrid and extract the 
               vector of data for the basin model for each layer.
               A bit odd because we only use the lookup function
               for the first point */
            lat=rad(latdeg);
            lon=rad(londeg);
            //elev /= 1000.0;
            double r=r0_ellipse(lat);
            Cartesian_point cp=basinmodel.gtoc(lat,lon,r);
            bool use_only_crust1p0(false);
            if(basinmodel.lookup(cp.x1,cp.x2,cp.x3))
            {
                cerr << "Lookup failed for sta "<<sta<<endl
                    << "Using Crust1.0 model for this station"<<endl;
                use_only_crust1p0=true;
            }
            int position[3];
            if(!use_only_crust1p0)
            {
            basinmodel.get_index(position);
            /* This is a one time program with a specialized known 
               model.  For that model and intended application a
               nearest neighbor approximation is adequate for 
               defining the basin model.  Well not really nearest
               neighbor but the lookup index position. */
            if(position[2]!=(basinmodel.n3-2))
            {
                cerr << "Lookup returned an unexpected value of "
                    << position[2] << " for station " <<sta<<endl
                    << "Expected to get "<<basinmodel.n3-2<<endl
                    << "Reset to that value - beware"<<endl;
                position[2]=basinmodel.n3-2;
            }
            /* now copy the basinmodel values to the output matrix*/
            int kk;
            //DEBUG
            /*
            cout << "Basin model only"<<endl;
            double zbsum;
            zbsum=basinmodel.depth(position[0],position[1],basinmodel.n3-1);
            zbsum+=5.0;   // top bad figure - debug only 
            */
            for(kk=0,k=basinmodel.n3-1;k>=0;++kk,--k)
            {
                for(int l=0;l<4;++l)
                    mergedmodel(kk,l)=basinmodel.val[position[0]][position[1]][k][l];
                //DEBUG
                /*
                if(kk!=0) zbsum+=mergedmodel(kk-1,3);
                cout << mergedmodel(kk,0) << " " 
                    << mergedmodel(kk,1) <<" "
                    << mergedmodel(kk,2) <<" "
                    << mergedmodel(kk,3)<<" "
                    <<zbsum<<endl;
                    */
            }
            //cout << endl;
            }  // End of if for accessing mergedmodel 
            k=basinmodel.n3-1;
            VelocityModel_1d Vp=c1p0.model(latdeg,londeg,string("Pvelocity"));
            VelocityModel_1d Vs=c1p0.model(latdeg,londeg,string("Svelocity"));
            VelocityModel_1d rho=c1p0.model(latdeg,londeg,string("Density"));
            //DEBUG
            //cout << "Crust1.0 model"<<endl;
            if(use_only_crust1p0)
            {
                cout << sta<<" "<<Vp.nlayers<<endl;
                for(int l=0;l<Vp.nlayers;++l)
                {
                    double dz;
                    if(l==(Vp.nlayers-1))
                    {
                        cout << Vp.v[l] << " "
                            << Vs.v[l] << " "
                            << rho.v[l] << " "
                            << 0.0<<endl;
                    }
                    else
                    {
                        dz=Vp.z[l+1]-Vp.z[l];
                        cout << Vp.v[l] << " "
                            << Vs.v[l] << " "
                            << rho.v[l] << " "
                            << dz <<endl;
                    }
                }
            }
            else
            {
            /* this positions kk to the upper crustal part of crust1.0.  
               Because we have to compute thickness not depths we have a 
               fairly elaborate calculation of thickness here. */
            kk=Vp.nlayers-4;
            /* reversed order of basinmodel makes this bottom of that model */
            double zbsed=basinmodel.depth(position[0],position[1],0);
            mergedmodel(k,0)=Vp.v[kk];
            mergedmodel(k,1)=Vs.v[kk];
            mergedmodel(k,2)=rho.v[kk];
            mergedmodel(k,3)=Vp.z[kk+1]-zbsed;
            //cout << "zbsed="<<zbsed<<endl<<endl;
            for(k=k+1,kk=kk+1;k<outsize;++k,++kk)
            {
                mergedmodel(k,0)=Vp.v[kk];
                mergedmodel(k,1)=Vs.v[kk];
                mergedmodel(k,2)=rho.v[kk];
                if(kk==(Vp.nlayers-1))
                    mergedmodel(k,3)=0.0;
                else
                    mergedmodel(k,3)=Vp.z[kk+1]-Vp.z[kk];
            }
            cout << sta << " " << outsize<<endl;
            cout << mergedmodel;
            }
            //DEBUG
            double zmoho;
            for(k=0,zmoho=0.0;k<outsize;++k)
                zmoho+=mergedmodel(k,3);
            cerr << "Crust10 moho depth="<<Vp.z[Vp.nlayers-1]<<" Sum of thickness="<<zmoho<<endl;
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
