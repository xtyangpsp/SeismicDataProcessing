/*
Modified from the original GLP code: testpwmigtt.
1. get vmodel file, delta, depth and other main variables from 
arguments.
2. purpose is to generate time-to-depth transformation function for 
the given model and source parameters.
3. write out to standard output (screen) and can be redirected to text file.
4. for the purpose of analyzing time-to-depth transform, this program freezes the
	ray-tracing type to be 'depth' or 'z'.
April 28, 2015
Xiaotao Yang
@Indiana University

2016.03.01, Xiaotao Yang: fixed bug, old version forgot to subtract scatter depth from EQUATORIAL_EARTH_RADIUS

*/
#include <string>
#include "ray1d.h"  //libseispp library
#include "tt.h"

/* //with option for ray tracing type.
void usage()
{
	cerr<<"getT2D vmodelfile distance depth [-rtt|--ray-tracing-type rttype]"<<endl
		<<"[ -rti|--ray-tracing-increment rtincrement][-md maxdepth][-mt maxtime]"<<endl
		<<"*Readin from mod1d format text file."<<endl
		<<"Defaults:"<<endl
		<<"    rttype: time; rtincrement: 0.1;"<<endl
		<<"    maxdepth: depth+5 km; maxtime: 200 seconds"<<endl
		<<"Units: distance: in great circle degrees; depth: in km;"<<endl
		<<"Options: -rtt: time or depth"<<endl;
	exit(-1);
}*/
//
void usage()
{
	cerr<<"getT2D vmodelfile distance depth "<<endl
		<<"[ -rti|--ray-tracing-increment rtincrement][-md maxdepth][-mt maxtime]"<<endl
		<<"*Readin from mod1d format text file."<<endl
		<<"Defaults:"<<endl
		<<"    rtincrement: 0.1 in depth;"<<endl
		<<"    maxdepth: depth+5 km; maxtime: 200 seconds"<<endl
		<<"Units: distance: in great circle degrees; depth: in km;"<<endl;
	exit(-1);
}
int main(int argc, char **argv)
{
    try{
    
    	if(argc<4) usage();
		string modname=string(argv[1]);
		double sdelta=atof(argv[2]);  //source distance in degree
		double sdepth=atof(argv[3]);  //source depth
		double maxdepth=sdepth+5.0;
		double maxtime(200.0);
		string rttype_tmp("depth"),rttype("z"); //raytracing type, either by time or by depth
		//default is time.
		double rtinc(0.1); //raytracing increment.

		int i;
		//nhdrs: number of header lines at the top.
		for(i=4;i<argc;++i)
		{
			string argstr(argv[i]);
			/*
			if(argstr=="-rtt" || argstr=="--ray-tracing-type")
			{
				++i;
				rttype_tmp=string(argv[i]);
			}
			else */
			if(argstr=="-rti" || argstr=="--ray-tracing-increment")
			{
				++i;
				rtinc=atof(argv[i]);
			}
			else if(argstr=="-md")
			{
				++i;
				maxdepth=atof(argv[i]);
			}
			else if(argstr=="-mt")
			{
				++i;
				maxtime=atof(argv[i]);
			}
			else
			{
				cerr<<"Illegal argument= "<<argstr<<endl;
				usage();
			}
		}
		if(rttype_tmp=="time") rttype="t";
		else if(rttype_tmp=="depth")
			rttype="z";
		else
		{
			cerr<<"Illegal rttype: "<<rttype_tmp<<endl;
			usage();
		}
        /* Frozen model  name for this test  ak135*/
        VelocityModel_1d Pvmod(modname,string("mod1d"),string("P"));
        VelocityModel_1d Svmod(modname,string("mod1d"),string("S"));
        /* set source delta (radians) and depth (km) for antelope calculator */
        //sdelta=sdelta*M_PI/180.0;  // 2 degrees
        double pslow,sslow;
        pslow=pphase_slowness(sdelta,sdepth);
        //sslow=sphase_slowness(sdelta,sdepth);
        RayPathSphere p_ray(Pvmod,pslow,maxdepth,maxtime,rtinc,rttype.c_str());
        RayPathSphere s_ray(Svmod,pslow,maxdepth,maxtime,rtinc,rttype.c_str());
        /* Now work down the computed ray path to compute comparisons */
        int np_p=p_ray.npts;
        //int itemp=np;
    
        //cout << "P ray comparison"<<endl
        //    << "delta(deg) depth t_iasp91 t_ray difference"<<endl;
        //cout<<"delta(deg)  depth(km)  Tp(seconds)  Ts(seconds)  Tps(seconds)"<<endl;
        cout<<"%delta(deg)  depth(km)  Tp(seconds)  Ts(seconds)  Tps(seconds)"<<endl;
        for(i=0;i<np_p;++i)
        {
        	double pzray,praydelta;
            double szray,sraydelta;
            pzray=p_ray.depth(i);
            szray=s_ray.depth(i);
            praydelta=p_ray.delta[i];
            sraydelta=s_ray.delta[i];
             // double time_s=s_ray.t[i]+(praydelta-sraydelta)*
//                           		pslow*EQUATORIAL_EARTH_RADIUS;
            double time_s=s_ray.t[i]+(praydelta-sraydelta)*
             		pslow*(EQUATORIAL_EARTH_RADIUS-fabs(pzray));
            double time_p=p_ray.t[i];
            double lag=time_s - time_p;
            		
        	
        	/*
        	double zray,raydelta;
        	zray=p_ray.depth(i);
            raydelta=p_ray.delta[i];
            pslow=pphase_slowness(raydelta,zray);
            sslow=sphase_slowness(raydelta,zray);
            RayPathSphere p_ray1(Pvmod,pslow,zray,maxtime,rtinc,rttype.c_str());
            RayPathSphere s_ray(Svmod,sslow,zray,maxtime,rtinc,rttype.c_str());
            int np_itmp=s_ray.npts;
            double time_s,depth_s, time_p,depth_p;
            
            if(i==0) 
            {time_s=s_ray.t[np_itmp-1];depth_s=s_ray.depth(np_itmp-1);}
            else 
            {time_s=s_ray.t[np_itmp-2];depth_s=s_ray.depth(np_itmp-2);}
            
            np_itmp=p_ray1.npts;
            if(i==0) 
            {time_p=p_ray1.t[np_itmp-1];depth_p=p_ray1.depth(np_itmp-1);}
            else 
            {time_p=p_ray1.t[np_itmp-2];depth_p=p_ray1.depth(np_itmp-2);}
            //debug
            //cout<<"pdepth: "<<zray<<", sdepth: "<<depth_s<<endl;
            //cout<<"ptime: "<<p_ray.t[i]<<", stime: "<<time_s<<endl;
            */
            printf("%10.5f  %8.3f  %10.5f  %10.5f  %10.5f\n",
            	praydelta*180.0/M_PI,pzray,time_p,time_s,lag);
        }
    }catch(SeisppError& serr)
    {
        cerr << serr.what()<<endl;
    }
}


