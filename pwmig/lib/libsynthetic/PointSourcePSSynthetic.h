#ifndef _POINTSOURCEPSYNTHETIC_H_
#define _POINTSOURCEPSYNTHETIC_H_
#include "SyntheticSeismogram.h"
#include <vector>
#include <VelocityModel_1d.h>
#include <ray1d.h>
#include "vectorcls.h"
#include "SphericalRayPathArray.h"

namespace SEISPP {
	using namespace std;
	using namespace SEISPP;

	class PointSourcePSSynthetic : public SyntheticSeismogram
	{
	public:
		PointSourcePSSynthetic(VelocityModel_1d& vsmods, VelocityModel_1d& vpmods, Pf *pf);
		//~PointSourcePSSynthetic();
		TimeSeries ComputeScalar(int nsamp, double dt, Hypocenter& hypo,
			double rlat, double rlon, double relev,string type);
		TimeSeries ComputeScalar(const TimeSeries& parent,
			Hypocenter& hypo,
			double rlat, double rlon, double relev,string type);
		ThreeComponentSeismogram Compute3C(int nsamp, double dt, Hypocenter& hypo,
			double rlat, double rlon, double relev,string units);
		ThreeComponentSeismogram Compute3C(const ThreeComponentSeismogram& parent,Hypocenter& hypo, double rlat, double rlon, double relev,string units);
		void P_scatter(double pointdepth, double plat, double hypodepth, double delta , double& Pinc_time, double& Pinc_theta);
		// Pinc_* are return values of P_scatter()
		double getlocalAZ(Hypocenter& hypo, Geographic_point& point, const double& rlat, const double& rlon);
		vectorcls getENZ(Hypocenter& hypo, Geographic_point& point, const double& rlat, const double& rlon, const double& relev, const double& Pinc_theta, const double& razimuth, const double& pertamp);
		int initPtime4event(Hypocenter& hypo);
		double ScatterPatternP2S(double costheta, double rhoperturb, double vsperturb, double vp, double vs);
	private:
		SphericalRayPathArray rays;
		/* purists might want this to be a single object, but here
		I use two parallel vectors to avoid this unnecessary
		complexity for a private variable.   Could easily change
		implementation of that proved essential */
		vector<Geographic_point> points;
		vector<double> amp, Ptimearray, Panglearray;
		//for Ptimearray and Panglearray, the initialization should be done
		//in the constructor.
		VelocityModel_1d vsmodel, vpmodel;
		
		//defines density perturbation (rhoperturb) and S velocity perturbation (vsperturb)
		double rhoperturb, vsperturb;
		int enableScatterPattern;
		///////////////////////////////////////////////////////////////////

		/* If a station is located beyond this t cutoff distance (radians)
		will not try to compute the point source response.  This is
		essential because we don't handle turning rays and for
		the for which I wrote this is links with the concept of 
		cutoff aperture in pwstack and (in progress) pwdecon*/
		double distance_cutoff;
		double vp0, vs0;// P and S velocity at free surface
		int rotationtype;
	};
	//now declare the functions in SEISPP namespace (contributed by xinliu)
	double triangularsolve(dmatrix& dm);
	vectorcls sphere2cardinal(const double& lat, const double& lon, const double& radius);
	double get_SSlowness(double degdelta, double depth);
	double phi_deltaS(const double& degdelta, const  double& ds, const double& sourcez, const double& pointz);
	dmatrix& vec2mat( vectorcls& vec, dmatrix& mat, int row0, int col0, int roworcol);
        dmatrix& vecAddTomat( vectorcls& vec, dmatrix& mat, int row0, int col0, int roworcol);
	vectorcls& mat2vec(dmatrix& mat, vectorcls& vec, int row0, int col0, int roworcol);
	//vectorcls mat2vec(dmatrix& mat, int row0, int col0, int roworcol);
}  // End SEISPP namespace encapsulation
#endif
