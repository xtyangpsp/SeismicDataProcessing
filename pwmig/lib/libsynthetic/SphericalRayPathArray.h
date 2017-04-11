#ifndef _SPHERICALRAYPATHARRAY_H_
#define _SPHERICALRAYPATHARRAY_H_
//#include "SyntheticSeismogram.h"
#include <vector>
#include <VelocityModel_1d.h>
#include <ray1d.h>
//#include "vectorcls.h"

namespace SEISPP {
using namespace std;
using namespace SEISPP;
class SphericalRayPathArray
{
    public:
        //SphericalRayPathArray(VelocityModel_1d& vmod,vector<double>& p,double zminin, double zmax, double dz);
	SphericalRayPathArray(VelocityModel_1d& vmod, vector<double>& rayp, double zminin, double zmax, double dzin);
        SphericalRayPathArray(VelocityModel_1d& vmod,double p0, double dp, int np, double zminin, double zmax, double dz);
        SphericalRayPathArray(const SphericalRayPathArray& parent);
        SphericalRayPathArray(){
		//blank constructor, do nothing here
	}
	
	pair<double,double> time_and_amp(double delta, double z);
      	SphericalRayPathArray& operator=(const SphericalRayPathArray& parent);
    private:
        int nrays;  // size of rays, cached for efficiency
        vector<RayPathSphere> rays;
        vector<double> slow;  // vector of ray parameters (s/km)
        /* zmin is the minimum depth allowed and a reference depth used for amplitude.
           That is, the amplitude for a ray from depth zmin to zero offset is given a unit
           amplitude factor */
        double zmin;  
        int nz;  // nominal depth steps that can be dependent on
        double dz;  // rays are made in equal depth steps, save here.
        //dmatrix amps;  // nz x nrays matrix of amplitudes derived from rays
	vector<vector<double> >  amps;// new amps based on ray parameter grid.
        double zfloor;  // actual depth floor = dz*(nz-1) cached for efficiency
};


}  // End SEISPP namespace encapsulation
#endif
