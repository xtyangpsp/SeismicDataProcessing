#include "coords.h"
#include "tt.h"
#include "interpolator1d.h"
#include "seispp.h"
#include "SimplePSPrimarySynthetic.h"
namespace SEISPP {
using namespace std;
using namespace SEISPP;
using namespace INTERPOLATOR1D;
SimplePSPrimarySynthetic::SimplePSPrimarySynthetic(Metadata& md,
        vector<double> scd, vector<double> sca)
{
    try{
        const string base_error("SimplePSPrimarySynethetic constructor:  ");
        sc_depth=scd;
        nlayers=sc_depth.size();
        if(sca.size()!=nlayers)
        {
            throw SeisppError(base_error
                    + "layer depth and amplitude input arrays have different sizes");
        }
        sc_amp=sca;
        vp0=md.get_double("surface_P_velocity");
        vs0=md.get_double("surface_S_velocity");
        ddelta=md.get_double("SimplePSPrimarySynthetic::ddelta");
        ndelta=md.get_int("SimplePSPrimarySynthetic::ndelta");
        if((ndelta*ddelta)>30.0) 
            throw SeisppError(base_error
                    + "ndelta * ddelta > 30 degrees not sensible - modify one or the other");
        tp=dmatrix(ndelta,nlayers);
        ts=dmatrix(ndelta,nlayers);
        pslow=dmatrix(ndelta,nlayers);
        sslow=dmatrix(ndelta,nlayers);
        int i,j;
        for(j=0;j<nlayers;++j)
        {
            double z=sc_depth[j]; 
            double delta;
            for(i=0,delta=0.0;i<ndelta;++i,delta+=ddelta)
            {
                pslow(i,j)=pphase_slowness(delta,z);
                sslow(i,j)=sphase_slowness(delta,z);
                tp(i,j)=pphasetime(delta,z);
                ts(i,j)=sphasetime(delta,z);
            }
        }
    }catch(...){throw;};
}
/* This private method gets the travel time by using the tables in
   the private area of this object */
double SimplePSPrimarySynthetic::delay_time(double deldeg, double sz, 
        int layer)
{
    try {
        double thispslow=pphase_slowness(deldeg, sz);
        int ip, is;
        ip=irregular_lookup(thispslow,pslow.get_address(0,layer),ndelta);
        is=irregular_lookup(thispslow,sslow.get_address(0,layer),ndelta);
        if((ip>=(ndelta-1)) || (is>=(ndelta-1)))
                throw SeisppError(
                    string("SimplePSPrimarySynethetic::delay_time")
                    + "  lookup overflowed internal arrays.  increase ndelta");

        const double deg2km(111.32);
        double ptdelta=thispslow*ddelta*static_cast<double>(ip-is)*deg2km;
        double dtp,dts;
        dtp=linear_scalar(thispslow,pslow(ip,layer),tp(ip,layer),
                pslow(ip+1,layer),tp(ip+1,layer));
        dts=linear_scalar(thispslow,sslow(is,layer),ts(is,layer),
                sslow(is+1,layer),ts(is+1,layer));
        double tau=dts+ptdelta-dtp;
        return(tau);
    }catch(...){throw;};
}
TimeSeries SimplePSPrimarySynthetic::ComputeScalar(const TimeSeries& parent,
        Hypocenter& hypo,double rlat, double rlon, double relev,string type)
{
    cout << "ComputeScalar procedure not implemented"<<endl;
    exit(-1);
}
TimeSeries SimplePSPrimarySynthetic::ComputeScalar(int nsamp,
        double dt, Hypocenter& hypo,
        double rlat, double rlon, double relev,string type)
{
    // This method can build a parent TimeSeries object, call the above, and
    // return it.  
    cout << "ComputeScalar plain procedure not implemented"<<endl;
    exit(-1);
}
ThreeComponentSeismogram SimplePSPrimarySynthetic::Compute3C(
        int nsamp, double dt,
        Hypocenter& hypo, double rlat, double rlon, double relev,
        string units)
{
    cout << "Compute3C plain procedure not implemented"<<endl;
    exit(-1);
}
ThreeComponentSeismogram SimplePSPrimarySynthetic::Compute3C(
        const ThreeComponentSeismogram& parent,Hypocenter& hypo,
        double rlat, double rlon, double relev,string units)
{
    const string base_error("SimplePSPrimarySynthetic::Compute3C method:  ");
    try {
        ThreeComponentSeismogram result(parent);
        result.u.zero();
        double delta,backaz,azimuth,phi;
        dist(rlat,rlon,hypo.lat,hypo.lon,&delta,&backaz);
	//dist returns backazimuth, convert to propagation azimuth
	azimuth=backaz+M_PI; 
	// phi is spherical coordinate angle from x1=E
	phi=M_PI_2-azimuth;
        /* shift to relative time if necessary */
        if(result.tref == absolute)
        {
            double patime=hypo.ptime(rlat,rlon,relev);
            result.ator(patime);
        }
        double deldeg=deg(delta);
        int i;
        for(i=0;i<nlayers;++i)
        {
            double tau=this->delay_time(deldeg,hypo.z,i);
            int it=result.sample_number(tau);
	    if(it>=0 && it<result.ns)
            	result.u(1,it)=sc_amp[i];
        }
        /* Now we have to set the transformation matrix for TRL coordinates.*/
        double slow0=pphase_slowness(deldeg,hypo.z);  // slightly inefficient to recompute this 
        SphericalCoordinate sc=PMHalfspaceModel(vp0,vs0,slow0*cos(phi), slow0*sin(phi));
        /* In a ThreeComponentSeismogram tmatrix is the matrix used as the
	forward transform to produce current data.  Data here are as if they
	were transformed to ray coordinates.  Hence, tmatrix code here 
	is identical to rotate method in ThreeComponentSeismogram object. */
        double a,b,c,d;
	a=cos(azimuth);
	b=sin(azimuth);
        c=cos(sc.theta);
        d=sin(sc.theta);
	result.tmatrix[0][0] = a;
	result.tmatrix[1][0] = b*c;
	result.tmatrix[2][0] = b*d;
	result.tmatrix[0][1] = -b;
	result.tmatrix[1][1] = a*c;
	result.tmatrix[2][1] = a*d;
	result.tmatrix[0][2] = 0.0;
	result.tmatrix[1][2] = -d;
	result.tmatrix[2][2] = c;
        // Perhaps unnecessary, but best to be sure these are set
        result.components_are_cardinal=false;
        result.components_are_orthogonal=true;
        result.rotate_to_standard();
        return(result);
    }catch(...){throw;};
}

}  // End SEISPP namespace encapsulation
