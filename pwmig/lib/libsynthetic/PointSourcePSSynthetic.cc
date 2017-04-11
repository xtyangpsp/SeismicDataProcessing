/* This object is the core of the implementation of this synthetic generator */
#include <iostream>
#include <fstream>
#include <cmath>
#include "coords.h"
#include "tt.h"
#include "stock.h"
#include "db.h"
#include "interpolator1d.h"
#include "seispp.h"
#include "PointSourcePSSynthetic.h" 
#include "vectorcls.h"
#include "SphericalRayPathArray.h"
//#include <Hypocenter.h>
namespace SEISPP {
	using namespace std;
	using namespace SEISPP;
	const double FLT_EPSILON=1e-19;

	//PointSourcePSSynthetic::PointSourcePSSynthetic(VelocityModel_1d& vmod,Pf *pf)

	PointSourcePSSynthetic::PointSourcePSSynthetic(VelocityModel_1d& vsmods, VelocityModel_1d& vpmods, Pf *pf)
	{
		const string base_error("PointSourcePSSynthetic constructor:  ");
		try{
			Tbl *t;
			t=pfget_tbl(pf,const_cast<char *>("ray_parameter_grid"));
			if(t==NULL) throw SeisppError(base_error
				+ "Missing required parameter ray_parameter_grid");
			int i;
			vector<double> raypvector;
			for(i=0;i<maxtbl(t);++i)
			{
				double rayp;
				char *line=(char *)gettbl(t,i);
				istringstream ss(line);
				ss >> rayp;
				raypvector.push_back(rayp);
			}
			freetbl(t,0);
			Metadata control(pf);
			double zmin=control.get_double("zmin");
			double zmax=control.get_double("zmax");
			double dz=control.get_double("dz");
			vp0=control.get_double("surface_P_velocity");
			vs0=control.get_double("surface_S_velocity");
			rotationtype=control.get_int("rotation_type");//1: curve; 0: line
			rays=SphericalRayPathArray(vsmods,raypvector,zmin,zmax,dz);
			vsmodel=vsmods;
			vpmodel=vpmods;
			//vmodel is a private member in this class
			//use the member function VelocityModel_1d::getv(double zin)
			distance_cutoff=control.get_double("distance_cutoff");
			//t=pfget_tbl(pf,"point_sources");
			//if(t==NULL) throw SeisppError(base_error + "Missing required parameter point_sources");
			string pointsrcf=control.get_string("pointsrcfile");
			int numpts=control.get_int("num_point_scatterers");
			points.reserve(numpts);
			amp.reserve(numpts);
			double r0;
			ifstream poisrcf(pointsrcf.c_str());
			//if((fp=fopen(pointsrcf.c_str(), "r"))){
			//	cout<<"The point source loc file doesn't exist!i\n";
			//}
			string line;
			while(getline(poisrcf, line))
			{
				//char *line=(char *)gettbl(t,i);
				double lat,lon,depth,r,ampin;
				Geographic_point gp;
				istringstream ss(line);
				ss >> lat;
				ss >> lon;
				ss >> depth;
				ss >> ampin;
				gp.lat=rad(lat);
				gp.lon=rad(lon);
				r0=r0_ellipse(lat);
				gp.r=r0-depth;
				points.push_back(gp);
				amp.push_back(ampin);


			}

			rhoperturb=control.get_double("density_perturb"); 
			vsperturb=control.get_double("Vs_perturb");
			// Notice: your migsimulation.pf file MUST contain density_perturb and Vs_perturb.
			enableScatterPattern=control.get_int("Apply_Scatter_Pattern");
			//fclose(poisrcf);
		} catch(...){throw;};
	}

	TimeSeries PointSourcePSSynthetic::ComputeScalar(int nsamp, double dt, Hypocenter& hypo,
		double rlat, double rlon, double rz,string type)
	{
		throw SeisppError(string("ComputeScalar method not implemented"));
	}
	TimeSeries PointSourcePSSynthetic::ComputeScalar(const TimeSeries& parent,
		Hypocenter& hypo,double rlat, double rlon, double rz,string type)
	{
		throw SeisppError(string("ComputeScalar method not implemented"));
	}
	ThreeComponentSeismogram PointSourcePSSynthetic::Compute3C(int nsamp, double dt, Hypocenter& hypo,
		double rlat, double rlon, double rz,string units)
	{
		throw SeisppError(string("Compute3C plain method not implemented"));
	}
	double triangularsolve(dmatrix& dm){

		int i,j,k;
		int nrows, ncols,tmp;//ncol is the number of columns
		int *piv;
		nrows=dm.rows();
		ncols=dm.columns();
		piv=new int[nrows];
		for(i=0;i<nrows;i++){
			piv[i]=i;
		}
		double* colmem=0;
		int maxind=0;
		int sign_cosntheta=1;// the sign of cos<n, theta>
		for(i=0;i<nrows-1;i++){
			colmem=dm.get_address(0,i);
			maxind=i;
			for(j=i+1;j<nrows;j++){
				//maxind=i; start from dm(i,i)
				if(colmem[piv[j]]>colmem[piv[maxind]]) 
					maxind=j;

				//	tmp=piv[i];piv[i]=piv[j];piv[j]=tmp;
				//swap(piv[i],piv[j]);
			}
			//swap (maxind+i).th and i.th row
			if(maxind!=i){
				tmp=piv[maxind];piv[maxind]=piv[i];piv[i]=tmp;
				sign_cosntheta= -sign_cosntheta;
				//when swap rows, determinant*(-1)
			}
			double ratio;
			for(j=i+1;j<nrows;j++){ //j is set to relative value of i.
				colmem=dm.get_address(0,i);
				ratio=-colmem[piv[j]]/colmem[piv[i]];
				//j-i is the j.th row, while 0 is the i.th row.
				for(k=0;k<ncols-i;k++){
					colmem[piv[j]+k*nrows]+=ratio*colmem[piv[i]+k*nrows];
				}
			}

		}
		//colmem=dm.get_address(piv[nrows-1],nrows-1);
		//double tannm=colmem[nrows]/colmem[0];//tan(Az-pi/2)
		colmem=dm.get_address(0,0);
		double detA=colmem[piv[0]];
		for(i=1;i<nrows;i++){
			detA*=colmem[piv[i]+i*nrows];

		}
		detA*=sign_cosntheta;
		colmem=dm.get_address(piv[nrows-1],nrows-1);
		delete[] piv;
		if(fabs(detA)<10e-8){
			if(colmem[nrows]>0)
				return M_PI/2;
			else
				return M_PI*1.5;

		}
		else{
			double deg1=0;	
			double tannm=colmem[nrows]/colmem[0];
			if(detA>0)
				deg1=atan(tannm);
			else
				deg1=atan(tannm)+M_PI;
			cout<<"This is detA: "<<detA<<endl;
			return deg1+M_PI/2;// this is the azimuth of S wave polarization.
		}
		//delete[] piv;
		//return colmem[nrows]/colmem[0];//return Xn

		/*if(fabs(colmem[0])<10e-8)

		return 0;//failed in calculating the 
		else{
		for(k=nrows-nrows;k<ncols-nrows;k++)
		colmem[nrows*k]/=colmem[0];
		//colmem[nrows+nrows]/=colmem[0];
		return 1;
		}
		*/
		//delete[] piv;
	}

	double PointSourcePSSynthetic::getlocalAZ(Hypocenter& hypo, Geographic_point& point, const double& rlat, const double& rlon){
		//az is the local azimuth of P-to-S converted wave ray path
		double az=0;
		dmatrix Aext=dmatrix(3,4);
		Aext(0,0)=cos(hypo.lat)*cos(hypo.lon);
		Aext(1,0)=cos(hypo.lat)*sin(hypo.lon);
		Aext(2,0)=sin(hypo.lat);
		Aext(0,1)=cos(point.lat)*cos(point.lon);
		Aext(1,1)=cos(point.lat)*sin(point.lon);
		Aext(2,1)=sin(point.lat);
		Aext(0,2)=-sin(rlat)*cos(rlon);
		Aext(1,2)=-sin(rlat)*sin(rlon);
		Aext(2,2)=cos(rlat);
		//Aext(1,5)=cos(rlat)*cos(rlon);
		//Aext(2,5)=cos(rlat)*sin(rlon);
		//Aext(3,5)=sin(rlat);
		Aext(0,3)= -sin(rlon);
		Aext(1,3)=cos(rlon);
		Aext(2,3)=0;

		//Matrix [A b] is constructed!
		double X3;
		//1=new double[3];
		az=triangularsolve(Aext);

		if(SEISPP_verbose) cout<<"The Azimuth of converted S polarization is:"<<deg(az)<<endl;
		//Aext(1,1)=sin()*cos();
		//az=X1[2];

		return az;
	}

	vectorcls sphere2cardinal(const double& lat, const double& lon, const double& radius){
		//lat==theta, lon==phi

		return vectorcls(radius*cos(lat)*cos(lon),radius*cos(lat)*sin(lon),radius*sin(lat));
	}
	// vec2mat() and mat2vec() are two functions implemented in SEISPP namespace
	// to convert between dmatrix and vectorcls objects.
	dmatrix& vec2mat( vectorcls& vec, dmatrix& mat, int row0, int col0, int roworcol){
		int i;
		double* mt=mat.get_address(row0, col0);
		double* vt=vec.get_address();
		if(roworcol){//col priority
			for(i=0;i<3;i++){
				mt[i]=vt[i];
			}

		}
		else{
			int nrows=mat.rows();
			for(i=0;i<3;i++){
				mt[i*nrows]=vt[i];
			}

		}

		return mat;

	}
        dmatrix& vecAddTomat( vectorcls& vec, dmatrix& mat, int row0, int col0, int roworcol){
                int i;
                double* mt=mat.get_address(row0, col0);
                double* vt=vec.get_address();
                if(roworcol){//col priority
                        for(i=0;i<3;i++){
                                mt[i]+=vt[i];
                        }

                }
                else{
                        int nrows=mat.rows();
                        for(i=0;i<3;i++){
                                mt[i*nrows]+=vt[i];
                        }

                }

                return mat;

        }

	vectorcls& mat2vec(dmatrix& mat, vectorcls& vec, int row0, int col0, int roworcol){
		// this function convert a row or a column of a dmatrix obj to 
		// a vectorcls obj.
		// either row0 or col0 should start from 0.
		// roworcol==1: col;
		// roworcol==0: row;
		int i;
		double* mt=mat.get_address(row0, col0);
		//vectorcls vec;
		double* vt=vec.get_address();
		if(roworcol){ // col priority
			//vectorcls vec(mt[0],mt[1],mt[2]);
			vt[0]=mt[0];vt[1]=mt[1];vt[2]=mt[2];
			return vec;

		}
		else{	
			int nrows=mat.rows();
			//vectorcls vec(mt[0],mt[nrows],mt[nrows*2]);
			vt[0]=mt[0];vt[1]=mt[nrows];vt[2]=mt[nrows*2];
			return vec;
		}
		//return vec;
	}

	void getlocalcardinalbase(dmatrix& mat, const double& lat, const double& lon){
		// produce the local cardinal coordinate system
		double* dm=mat.get_address(0,0);
		dm[0]= -sin(lon);
		dm[1]= cos(lon);
		dm[2]= 0;
		dm[3]= -sin(lat)*cos(lon);
		dm[4]= -sin(lat)*sin(lon);
		dm[5]= cos(lat);
		dm[6]= cos(lat)*cos(lon);
		dm[7]= cos(lat)*sin(lon);
		dm[8]= sin(lat);



	}
	double get_SSlowness(double degdelta, double depth){
		//notice that degdelta in in degrees, depth in km.
		extern char    *optarg;
		extern int      optind;
		void           *library ;
		int             c, errflg = 0;
		const int maxsize=8;
		char method[maxsize]="tttaup\0";
		char model[maxsize]="iasp91\0";
		char phase[maxsize]="S\0";
		TTGeometry geometry;
		Tbl *slow=0;
		Hook *hook=0;
		int mode=2;
		geometry.source.lat=0.0;
		geometry.source.lon=0.0;
		geometry.source.z=depth;
		geometry.source.time = 0.0 ;
		strcpy(geometry.source.name, "SOURCE" ) ;

		geometry.receiver.lat=0.0;
		geometry.receiver.lon=0+degdelta;
		geometry.receiver.z=0.0;
		geometry.receiver.time = 0.0 ;
		strcpy(geometry.receiver.name, "RECEIVER" ) ;

		//now all initial values & parameters are defined.
		int result=ucalc(method, model, phase, mode, &geometry, &slow, &hook);

		if(result==-1 && SEISPP_verbose) cout<<"result=-1, no solutions found\n";
		TTSlow *u;
		u=(TTSlow *) gettbl(slow, 0);
		//cout<<"phase="<<u->phase<<endl;

		double u_horizontal=sqrt((u->ux*u->ux)+(u->uy)*(u->uy));
		//cout<<"S wave slowness="<<u_horizontal<<endl;
		freetbl(slow,0);
		//delete u;
		free_hook(&hook);
		//freetbl(slow,0);
		//free(u);
		return  u_horizontal;
	}

	vectorcls PointSourcePSSynthetic::getENZ(Hypocenter& hypo, Geographic_point& point, const double& rlat, const double& rlon, const double& relev, const double& Pinc_theta, const double& razimuth, const double& pertamp){
		// basically, this function calculate the direction of S wave polarization
		// in the local geographic coordinate system at receiver;
		// the amplitude of S polar is normalized to 1
		double pazimuth, delta;
		//dist(rlat,rlon,points[i].lat,points[i].lon,&delta,&azimuth);
		dist(point.lat,point.lon,rlat,rlon,&delta,&pazimuth);//actually, this is the azimuth of converted S wave @ the point scatterer
		//reverse source & receiver to get the azimuth @ the point source.
		double deldeg=deg(delta);
		double pointdepth=r0_ellipse(deg(point.lat))-point.r;//r0_ellipse(degree)
		//double slow0=pphase_slowness(deldeg,pointdepth);//
		double slow0=get_SSlowness(deldeg,pointdepth);
		if(SEISPP_verbose) cout<<"slow0="<< slow0<<endl;
		//double slow0=sphase_slowness(10,pointdepth/10);

		//P phase slowness or S phase slowness?? assume P polarization is perpendicular to S?
		/*double slow0=pphase_slowness(deldeg,r0_ellipse(point.lat)-point.r);
		SphericalCoordinate sc=PMHalfspaceModel(vp0,vs0,slow0*sin(azimuth), slow0*cos(azimuth));  */

		double vs=vsmodel.getv(pointdepth);//get S velocity from 1D velocitymodel
		double vp=vpmodel.getv(pointdepth);
		double theta,theta1;
		if(slow0==-1 || deldeg==0){//the distance (deldeg) is tooooo small, treat the path as vertical
			theta=theta1=0;
		}
		else{
			theta=asin(slow0*vs0);//slow0/(1/v0)
			theta1=asin(r0_ellipse(deg(rlat))*vs*sin(theta)/(vs0*point.r));//S scattered
			// theta1 is the elevation angle of scattered S wave at depth of the point
			if(SEISPP_verbose) cout<<"theta of S wave @ surface is "<< deg(theta)<<endl;
			if(SEISPP_verbose)  cout<<"while delta=  degree " <<deldeg<<endl;
		}
		/*SlownessVector pslowvector=hypo.pslow(point.lat, point.lon, pointdepth);
		double azPin=pslowvector.azimuth();
		double thetap=asin(pslowvector.mag()*vp);*/
		double azPin;
		dist(point.lat,point.lon,hypo.lat,hypo.lon,&delta,&azPin);
		azPin-=M_PI;
		vectorcls Pin=sphere2cardinal(M_PI/2-Pinc_theta, M_PI/2-azPin, 1);
		//Pin expressed in local geographical coordinate system
		vectorcls Lp=sphere2cardinal(M_PI/2-theta1,M_PI/2-pazimuth,1);
		//expressed in local geographical coordinate system
		dvector dv(3);//3 component vector derived from dmatrix class
		vec2mat(Lp, dv, 0, 0, 1);
		//first two 0s indicate the starting row and col number
		//the last parameter 1 indicates column major
		dmatrix Atrans(3,3);
		getlocalcardinalbase(Atrans, point.lat, point.lon);
		dv=Atrans*dv;
		mat2vec(dv, Lp, 0, 0, 1);
		//-----------------rotate Pin to World coordinate system---------
		vec2mat(Pin, dv, 0, 0, 1);
		dv=Atrans*dv;
		//the same transform matrix Atrans
		mat2vec(dv, Pin, 0, 0, 1);
		//-----------------get Rsc from Pin and Lp----------------
		vectorcls Rsc=Lp*(Pin*Lp); // Rsc: S wave polarization@point
		Rsc/=modul(Rsc);

		if(deldeg==0){
			if(SEISPP_verbose) cout<<"hahah"<<endl;
			double cos_sc_angle=dotproduct(Pin, Lp);
			vectorcls spolar= sphere2cardinal(0, M_PI/2-azPin, 1);
			if(enableScatterPattern) spolar*=ScatterPatternP2S(cos_sc_angle, rhoperturb, vsperturb, vp, vs);
			return spolar;
			
		}
		//vectorcls focal=sphere2cardinal(hypo.lat, hypo.lon, 1);
		//vectorcls n1=sphere2cardinal(hypo.lat, hypo.lon, 1)*sphere2cardinal(point.lat, point.lon,1);// this is a cross product of 2 vectors
		//n1=fXp/|fXp|
		//n1/=modul(n1);//scale to a modulus of one
		//vectorcls Rsc=n1*Lp;//also cross product, Lp==Lsc
		//Rsc/=modul(Rsc);//scale to a modulus of one
		vectorcls Tp=sphere2cardinal(point.lat, point.lon,1)*Lp;
		Tp/=modul(Tp);
		vectorcls Rp=Lp*Tp;
		Rp/=modul(Rp);
		double* ptr=dv.get_address(0,0);
		ptr[0]=dotproduct(Rsc, Tp);// minus because the transform is in Left handed 
		ptr[1]=dotproduct(Rsc, Rp);// coords. now I drop this minus sign...
		ptr[2]=0;

		getlocalcardinalbase(Atrans, M_PI/2-theta, M_PI/2-razimuth);
		/*
		dmatrix Btrans(3,3);//local cardinal coordinate system @ receiver
		getlocalcardinalbase(Btrans, rlat, rlon);
		Btrans=tr(Btrans);
		//here you should define Rpnew and Lpnew, which 
		//correspond to the ray path bases at the receiver.
		vec2mat(Tp, Atrans, 0, 0, 1);
		vec2mat(Rpnew, Atrans, 0, 1, 1);
		vec2mat(Lpnew, Atrans, 0, 2, 1);

		dv=Btrans*Atrans*dv;
		*/
		dv=Atrans*dv;
		vectorcls Spolar;
		mat2vec(dv, Spolar, 0, 0, 1);
		double cos_sc_angle=dotproduct(Pin, Lp);
		if(enableScatterPattern) Spolar*= ScatterPatternP2S(cos_sc_angle, rhoperturb, pertamp*vsperturb, vp, vs); //scalar product
		//double scat_pat=ScatterPatternP2S(0.1, rhoperturb, vsperturb, vp, vs);
		//cout<<"scatter_pattern="<<scat_pat<<endl;
                //if(enableScatterPattern) Spolar*=0.0018;
		//ScatterPatternP2S(0.1 , rhoperturb, vsperturb, vp, vs); //scalar product

		//above line scales the polarization by scattering pattern.
		//cout<<"scatter_angle="<<deg(acos(cos_sc_angle))<<", Spolar="<<Spolar<<endl;
		return Spolar;
		//return vectorcls(1,2,3);//	
	}

	double getPphasetime(double a,double b)
	{
		TTGeometry      geometry ;
		geometry.source.lat = 0.0 ;
		geometry.source.lon = 0.0 ;
		geometry.source.z = 0.0 ;
		geometry.source.time = 0.0 ;
		strcpy(geometry.source.name, "SOURCE" ) ;

		geometry.receiver.lat = 0.0 ;
		geometry.receiver.lon = 0.0 ;
		geometry.receiver.z = 0.0 ;
		geometry.receiver.time = 0.0 ;
		strcpy(geometry.receiver.name, "RECEIVER" ) ;
		geometry.source.z = b ;
		geometry.receiver.lon = a ;
		Tbl            *times1 = 0;
		Hook           *hook = 0 ;
		ttcalc ( const_cast<char *>("tttaup"), const_cast<char *>("iasp91"), const_cast<char *>("P"), 2, &geometry, &times1, &hook ) ;
		TTTime *atime ;
		atime = (TTTime *) gettbl(times1, 0);
		return atime->value;
	}
	double getPnphasetime(double a,double b)
	{
		TTGeometry      geometry ;
		geometry.source.lat = 0.0 ;
		geometry.source.lon = 0.0 ;
		geometry.source.z = 0.0 ;
		geometry.source.time = 0.0 ;
		strcpy(geometry.source.name, "SOURCE" ) ;

		geometry.receiver.lat = 0.0 ;
		geometry.receiver.lon = 0.0 ;
		geometry.receiver.z = 0.0 ;
		geometry.receiver.time = 0.0 ;
		strcpy(geometry.receiver.name, "RECEIVER" ) ;
		geometry.source.z = b ;
		geometry.receiver.lon = a ;
		Tbl            *times1 = 0;
		Hook           *hook = 0 ;
		ttcalc ( const_cast<char *>("tttaup"), const_cast<char *>("iasp91"), const_cast<char *>("Pn"), 2, &geometry, &times1, &hook ) ;
		TTTime *atime ;
		atime = (TTTime *) gettbl(times1, 0);
		return atime->value;
	}
	double getPphaseslowness(double a,double b)
	{
		TTGeometry      geometry ;
		geometry.source.lat = 0.0 ;
		geometry.source.lon = 0.0 ;
		geometry.source.z = 0.0 ;
		geometry.source.time = 0.0 ;
		strcpy(geometry.source.name, "SOURCE" ) ;

		geometry.receiver.lat = 0.0 ;
		geometry.receiver.lon = 0.0 ;
		geometry.receiver.z = 0.0 ;
		geometry.receiver.time = 0.0 ;
		strcpy(geometry.receiver.name, "RECEIVER" ) ;
		geometry.source.z = b ;
		geometry.receiver.lon = a ;
		Tbl            *slow1 = 0;
		Hook           *hook = 0 ;
		ucalc ( const_cast<char *>("tttaup"), const_cast<char *>("iasp91"), const_cast<char *>("P"), 2, &geometry, &slow1, &hook);
		TTSlow *u ;
		u = (TTSlow *) gettbl(slow1, 0) ;
		return u->ux;
	}

	double phi_deltaS(const double& degdelta, const  double& ds, const double& sourcez, const double& pointz){
	// note that delta and ds should be in degree
		
		return getPphaseslowness(degdelta+ds, sourcez)-getPphaseslowness(ds, pointz);
	} 

	void PointSourcePSSynthetic::P_scatter(double pointdepth, double plat, double hypodepth, double delta , double& Pinc_time, double& Pinc_theta){
		//return the last two parameters
		// basically, this function compute the travel time and incident angle
		// for teleseismic P wave at the point scatterer. All angles are in radians.
		int i;// iteration variable
		int maxn=22;// number of maximum iterations
		const double converge_eps=1e-8, epsi=1e-5;//previously, converge_eps was set to 1e-9, which is too small.
		delta=deg(delta);//this delta is a copy of the one in parent funciton
		//delta is converted to degree
		/*if(delta<30 || delta>99){
		cout<<"no converted phases! degdelta="<<delta<<", depth="<<pointdepth <<endl;
		Pinc_time=1e6;
		Pinc_theta=0;

		}*/

		//set up the initial value of deltaS (ds)
		double r0=r0_ellipse(deg(plat));
		//double dipping=(20-12)/(30-100)*(delta-30)+20;//
		double dipping=(10-1.0)/(30.0-100.0)*(delta-30)+10;
		dipping=rad(dipping);
		const double x0=deg(pointdepth*(dipping)/r0_ellipse(0));
		//const double x0=deg(pointdepth*(M_PI)/9/2/r0_ellipse(0));
		//assume the incident angle is not more than 20 degrees
		double xk=x0;
		double phi_xk;
		for(i=0;i<maxn;i++){
			phi_xk=phi_deltaS(delta,xk,hypodepth,pointdepth);
			xk-=phi_xk*epsi/(phi_deltaS(delta,xk+epsi,hypodepth,pointdepth)-phi_xk);
			if(fabs(phi_xk)<converge_eps) break;
		}
		if(i==maxn) {
			if(SEISPP_verbose) cout<<"Not converged! degdelta="
				<<delta<<", xk="<<xk 
				<<",phi_xk"<<phi_xk<<", dipping"<<deg(dipping)<<endl;

			Pinc_time=1e6;
			Pinc_theta=0;
		}
		else{
			//now we've got xk (in degree), or the approximation of deltaS
			double t1=getPphasetime(delta+xk,hypodepth);
			if(t1==-1)
			{
				t1=getPnphasetime(delta+xk,hypodepth);
			}
			if(t1==-1)
			{	
				cout<<"Travel Time Error!!"<<endl;
				exit(-1);
			}
			if(SEISPP_verbose) cout<<"t1 is computed: "<<t1<<endl;
			double t2=getPphasetime(xk,pointdepth);
			if(t2==-1) t2=getPnphasetime(xk,pointdepth);
			if(t2==-1)
			{	
				cout<<"Travel Time Error!!"<<endl;
				exit(-1);
			}
			if(SEISPP_verbose) cout<<"t2 is computed: "<<t2<<endl;
			Pinc_time=t1-t2;
			if(SEISPP_verbose) cout<<"Pinc_time is computed: "<<Pinc_time<<endl;
			double vp=vpmodel.getv(pointdepth);
			//double r0=r0_ellipse(plat);
			double sinetheta=getPphaseslowness(delta+xk,hypodepth)*vp0;
			Pinc_theta=asin(r0*vp*sinetheta/(vp0*(r0-pointdepth)));
		}
		//incorrect now
	}

	int PointSourcePSSynthetic::initPtime4event(Hypocenter& hypo){
		//initialize P time and angle arrays for a specific event.
		Ptimearray.clear(); Panglearray.clear();
		int npoints=points.size();
		for(int i=0;i<npoints;i++){
			double p_2scatter_time;
			double pointidepth=r0_ellipse(deg(points[i].lat))-points[i].r;
			double azPin,deltap, Pinc_theta;
			dist(points[i].lat,points[i].lon,hypo.lat,hypo.lon,&deltap,&azPin);
			azPin-=M_PI;
			if(deltap<rad(30)){
				if(SEISPP_verbose)  cout<<"no converted phases! degdelta="<<deg(deltap)<<", depth="<<pointidepth <<endl;
				p_2scatter_time=1e6;
				Pinc_theta=0;
			}
			else{		
				P_scatter(pointidepth, points[i].lat, hypo.z, deltap , p_2scatter_time, Pinc_theta);
				if(SEISPP_verbose) cout<<"P to scatterer time:"<<p_2scatter_time<<", point i:"
					<<i<<endl;
			}
			Ptimearray.push_back(p_2scatter_time);
			Panglearray.push_back(Pinc_theta);

		}
		return 0;

	}
	double PointSourcePSSynthetic::ScatterPatternP2S(double costheta, double rhoperturb, double vsperturb, double vp, double vs){
		// ACCORDING TO OUR DEFINITION, theta is the scattering angle
		// which is between [0, pi].
		double sintheta=sqrt(1-costheta*costheta);
		double sin2theta=2*sintheta*costheta;
		return -((vsperturb+rhoperturb)*(vs/vp)*sin2theta+rhoperturb*sintheta);
		//	return 0;
	}
	ThreeComponentSeismogram PointSourcePSSynthetic::Compute3C(const ThreeComponentSeismogram& parent, Hypocenter& hypo,double rlat, double rlon, double relev,string units)
	{
		//before calling this function, you should initialize the Ptimearray
		//and Panglearray vectors by calling PointSourcePSSynthetic::initPtime4event().
		try{

			ThreeComponentSeismogram result(parent);
			result.u.zero();
			double delta,azimuth;
			//dist(rlat,rlon,hypo.lat,hypo.lon,&delta,&azimuth);
			int npoints=points.size(); //get the size of vector points
			//the following line gets the P arrival time at receiver
			double patime=hypo.ptime(rlat,rlon,relev);
			if(result.tref == absolute)
			{
				double patime=hypo.ptime(rlat,rlon,relev);
				result.ator(patime);
			}


			//this loop calculate the amplitude & relative time for 
			// each scattering point
			for(int i=0;i<npoints;i++){
				//hook=0;
				double  s_time;
				double tau,amplitude;
				double pointidepth=r0_ellipse(deg(points[i].lat))-points[i].r;
				//using this dist() to get delta between scattering point and receiver
				dist(rlat,rlon,points[i].lat,points[i].lon,&delta,&azimuth);
				azimuth-=M_PI;//convert backazimuth to azimuth @ receiver
				// the return values are put into a std::pair struct;
				// note that rays is a SphericalRayPathArray object;
				//convert elevation of points[i] to depth.
				try{
					pair<double,double> tnamp=rays.time_and_amp(delta, pointidepth); //geometrical spreading term
					s_time=tnamp.first;//delta is in radians?
					if(SEISPP_verbose) cout<<s_time<<" "<<sphasetime(deg(delta), pointidepth) <<" "<<endl;
					if(delta*6370.0>distance_cutoff) continue;
					//if(s_time==0) throw SeisppError("distance exceeds cut off limit");
					if(s_time==0) continue;
					tau=Ptimearray[i]+s_time-patime;
					amplitude=tnamp.second;
					int it=result.sample_number(tau);
					if(it>=0 && it<result.ns){
						if(rotationtype){
							vectorcls spolar=getENZ(hypo, points[i], rlat, rlon, relev, Panglearray[i], azimuth, amp[i]);//the converted S wave azimuth@receiver
							spolar*=amplitude;
							//get_SSlowness(10 , pointidepth);
							//cout<<"S slow:"<<sphase_slowness(10, pointidepth)<<endl;
							//vec2mat(spolar, result.u, 0, it, 1);
                                                        vecAddTomat(spolar, result.u, 0, it, 1);
							//vec2mat effectively copies spolar to the 
							//it.th column of matrix result.u 		
						}   	 
						else{

							double a,b;
							double PolarAz=getlocalAZ(hypo,points[i],rlat,rlon);
							a=cos(PolarAz);
							b=sin(PolarAz);
							result.u(0,it)=amplitude*b;
							result.u(1,it)=amplitude*a;
						}

					}	//result.u(1,it)+=amplitude;//I don't know if it works this way.
				}
				catch(SeisppError& serr){
					serr.log_error();
					if(SEISPP_verbose) cout<<"an error just occurred!"<<endl;
					continue;
				}						    // It shouble be superposition of each waveform
				// scattered from a single point.

			}
			//	dist(rlat,rlon,hypo.lat,hypo.lon,&delta,&azimuth);
			//	double deldeg=deg(delta);//get delta between receiver and hypocenter in degree
			/* Now we have to set the transformation matrix for TRL coordinates.*/
			//        double slow0=pphase_slowness(deldeg,hypo.z);  // should deldeg be in a unit of degree?
			//        SphericalCoordinate sc=PMHalfspaceModel(vp0,vs0,slow0*sin(azimuth), slow0*cos(azimuth));
			/* This builds the transpose of the tmatrix used in the rotate procedure method in
			ThreeComponentSeismogram using ray coordinates.  Since the coordinate system is
			orthogonal this is effectively setting the inverse */
			/*       double a,b,c,d;
			a=cos(sc.phi);
			b=sin(sc.phi);
			c=cos(sc.theta);
			d=sin(sc.theta);
			result.tmatrix[0][0] = a*c;
			result.tmatrix[0][1] = -b;
			result.tmatrix[0][2] = a*d;
			result.tmatrix[1][0] = b*c;
			result.tmatrix[1][1] = a;
			result.tmatrix[1][2] = b*d;
			result.tmatrix[2][0] = -d;
			result.tmatrix[2][1] = 0.0;
			result.tmatrix[2][2] = c;
			// Perhaps unnecessary, but best to be sure these are set
			*/
			/*	double a,b;
			double PolarAz=getlocalAZ(hypo,points[i],rlat,rlon);
			a=cos(PolarAz);
			b=sin(PolarAz);
			//	result.tmatrix.zero();
			result.tmatrix[0][0] = a;
			result.tmatrix[0][1] = -b;
			result.tmatrix[0][2] = 0;
			result.tmatrix[1][0] = b;
			result.tmatrix[1][1] = a;
			result.tmatrix[1][2] = 0;
			result.tmatrix[2][0] = 0;
			result.tmatrix[2][1] = 0;
			result.tmatrix[2][2] = 0;
			*/
			result.components_are_cardinal=true;//already rotated!--Xin
			result.components_are_orthogonal=true;
			//      result.rotate_to_standard();
			return(result);


			/* shift to relative time if necessary */
			/*if(result.tref == absolute)
			{
			double patime=hypo.ptime(rlat,rlon,relev);
			result.ator(patime);
			} */
			/*for(i=0;i<nlayers;++i)
			{
			double tau=this->delay_time(deldeg,hypo.z,i);
			int it=result.sample_number(tau);
			if(it>=0 && it<result.ns)
			result.u(1,it)=sc_amp[i];
			}*/

		} catch(...){throw;};
	}

}//end of namespace SEISPP
