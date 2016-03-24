#include <string>
#include <fstream>
#include "perf.h"
#include "gclgrid.h"
#include "seispp.h"
using namespace std;
using namespace SEISPP;
void usage()
{
    cerr << "merge_vmodels model1 model2 outmodel"<<endl
        << " ** model1 is dominant. model2 will be replaced with model1 where model1 has coverage."<<endl;
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
    if(argc!=4) usage();
    try {
        int k,kk;
        //Crust1_0 c1p0;
        string modelfile1(argv[1]);
        string modelfile2(argv[2]);
        string modelfileout(argv[3]);
        GCLvectorfield3d vmodelin1(modelfile1);
        GCLvectorfield3d vmodelin2(modelfile2);
        /* This procedure is super bubba, but it is necessary because 
           the upper layers in this model are very very thin and the 
           lookup method used in the gclgrid library is not appropriate
           for this geometry.   */
        pad_top_model(vmodelin1);
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
            Cartesian_point cp=vmodel.gtoc(lat,lon,r);
            //bool use_only_crust1p0(false);
            if(vmodel.lookup(cp.x1,cp.x2,cp.x3))
            {
                cerr << "Lookup failed for sta "<<sta<<endl
                    << "Skipped this station"<<endl;
                //use_only_crust1p0=true;
                continue;
            }
            int position[3];
            //if(!use_only_crust1p0)
            //{
            vmodel.get_index(position);
            /* This is a one time program with a specialized known 
               model.  For that model and intended application a
               nearest neighbor approximation is adequate for 
               defining the basin model.  Well not really nearest
               neighbor but the lookup index position. */
            if(position[2]!=(vmodel.n3-2))
            {
                cerr << "Lookup returned an unexpected value of "
                    << position[2] << " for station " <<sta<<endl
                    << "Expected to get "<<vmodel.n3-2<<endl
                    << "Reset to that value - beware"<<endl;
                position[2]=vmodel.n3-2;
            }
            /* now copy the vmodel values to the output matrix*/
            int kk;
            //DEBUG
            /*
            cout << "Basin model only"<<endl;
            double zbsum;
            zbsum=vmodel.depth(position[0],position[1],vmodel.n3-1);
            zbsum+=5.0;   // top bad figure - debug only 
            */
            for(kk=0,k=vmodel.n3-1;k>=0;++kk,--k)
            {
                for(int l=0;l<4;++l)
                    modeldata(kk,l)=vmodel.val[position[0]][position[1]][k][l];
                //DEBUG
                /*
                if(kk!=0) zbsum+=modeldata(kk-1,3);
                cout << modeldata(kk,0) << " " 
                    << modeldata(kk,1) <<" "
                    << modeldata(kk,2) <<" "
                    << modeldata(kk,3)<<" "
                    <<zbsum<<endl;
                    */
            }
            //cout << endl;
            //}  // End of if for accessing modeldata 
            //k=vmodel.n3-1;
            /*
            VelocityModel_1d Vp=c1p0.model(latdeg,londeg,string("Pvelocity"));
            VelocityModel_1d Vs=c1p0.model(latdeg,londeg,string("Svelocity"));
            VelocityModel_1d rho=c1p0.model(latdeg,londeg,string("Density"));
            */
            //DEBUG
            //cout << "Crust1.0 model"<<endl;
            /* this positions kk to the upper crustal part of crust1.0.  
               Because we have to compute thickness not depths we have a 
               fairly elaborate calculation of thickness here. */
            //kk=Vp.nlayers-4;
            /* reversed order of vmodel makes this bottom of that model */
            /*double zbsed=vmodel.depth(position[0],position[1],0);
            modeldata(k,0)=Vp.v[kk];
            modeldata(k,1)=Vs.v[kk];
            modeldata(k,2)=rho.v[kk];
            modeldata(k,3)=Vp.z[kk+1]-zbsed;
            //cout << "zbsed="<<zbsed<<endl<<endl;
            for(k=k+1,kk=kk+1;k<outsize;++k,++kk)
            {
                modeldata(k,0)=Vp.v[kk];
                modeldata(k,1)=Vs.v[kk];
                modeldata(k,2)=rho.v[kk];
                if(kk==(Vp.nlayers-1))
                    modeldata(k,3)=0.0;
                else
                    modeldata(k,3)=Vp.z[kk+1]-Vp.z[kk];
            }
            */
            cout << sta << " " << outsize<<endl;
            cout << modeldata;

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
