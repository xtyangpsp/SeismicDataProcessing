#ifndef _SIMPLEPSPRIMARYSYNTHETIC_H_
#define _SIMPLEPSPRIMARYSYNTHETIC_H_
#include <vector>
#include "SyntheticSeismogram.h"

namespace SEISPP {
using namespace std;
using namespace SEISPP;
class SimplePSPrimarySynthetic : public SyntheticSeismogram
{
    public:
        SimplePSPrimarySynthetic(Metadata& md, vector<double> scd,
                vector<double> sca);
        TimeSeries ComputeScalar(int nsamp, double dt, Hypocenter& hypo,
                            double rlat, double rlon, double relev,string type);
        TimeSeries ComputeScalar(const TimeSeries& parent,
             Hypocenter& hypo,double rlat, double rlon, double relev,string type);
        ThreeComponentSeismogram Compute3C(int nsamp, double dt, Hypocenter& hypo,
                 double rlat, double rlon, double relev,string units);
        ThreeComponentSeismogram Compute3C(const ThreeComponentSeismogram& parent,
             Hypocenter& hypo, double rlat, double rlon, double relev,string units);
    private:
        /*size variables cached here for efficiency.  
         nlayers=number of impulses created
         ndelta = size of internal travel time curves for each layer*/
        int nlayers, ndelta;
        /* The internal travel times run from 0 to ndelta*ddelta.
         ddelta input is in km but internally is in degrees to allow
         it to work cleaner with the antelope travel time functions.*/
        double ddelta; 
        /* This algorithm is dumb and uses a fixed amplitude for all
           ray parameters.  i.e. it does not use physical conversion
           coefficients.  The amplitude (normalized to 1 for incident wave)
           is stored in sc_amp and sc_depth is a parallel vector storing
           depths.   nlayers is computed as size() of these vectors*/
        vector<double> sc_depth,sc_amp;
        /* These are ndelta by nlayers arrays to hold travel time for
           P, travel time for S, slowness P, and slowness S respectively */
        dmatrix tp,ts,pslow,sslow;
        double delay_time(double deldeg, double sz, int layer);
        double vp0,vs0;  // P and S velocity at free surface (km/s)
};
}  // End SEISPP namespace encapsulation
#endif
