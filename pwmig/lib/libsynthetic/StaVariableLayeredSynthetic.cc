#include <math.h>
#include <fstream>
#include "StaVariableLayeredSynthetic.h"
using namespace std;
using namespace SEISPP;
/* this is the fortran subroutine interface to the Kennett plane wave synthethic
   algorithm.*/
extern "C" void kntsyn_(int *nlyr, float *alfm, float *betm, float *rhom,
      float *thikm, float *pr, int *full,float *dt, int *numpts,
       float *tsigma, float *wlevel, float *u0, float *w0, float *u1,
         float *w1,float *tn, float *rfr, int *ierr);

/* Rather primitive function extracts the leaf name of a unix file
   and returns it.  For this code when driven by a list of file names
   this is assumed to be a seismic station name. */
string extract_sta_name(string path)
{
    string sta;
    size_t right;
    right=path.rfind("/");
    if(right==string::npos)
        sta=path;
    else
        sta.assign(path,right+1,string::npos);
    return(sta);
}
StaVariableLayeredSynthetic::StaVariableLayeredSynthetic(list<string> filelist,
        double tsig,double wlev,bool MakeRF)
{
    const string base_error("StaVariableLayeredSynthetic constructor:  ");
    tsigma=static_cast<float>(tsig);
    wlevel=static_cast<float>(wlev);
    ConvertToRF=MakeRF;
    try {
        LayeredModel thismod;
        list<string>::iterator mptr;
        for(mptr=filelist.begin();mptr!=filelist.end();++mptr)
        {
            string sta=extract_sta_name(*mptr);
            ifstream din;
            din.open(mptr->c_str(),ios::in);
            if(din.fail())
                throw SeisppError(base_error + "Open failed on model file "
                        + *mptr);
            char inputline[256];
            double Pvin,Svin,density,dzin;
            thismod.name=sta;
            while(din.getline(inputline,256))
            {
                stringstream ss(inputline);
                ss >> Pvin; ss>>Svin;  ss>>density; ss>>dzin;
                if(Pvin<=0.0 || Svin <= 0.0 || density<=0.0)
                    throw SeisppError(base_error
                            + "Illegal model parameters for mode file"
                            + *mptr);
                thismod.alpha.push_back(Pvin);
                thismod.beta.push_back(Svin);
                thismod.rho.push_back(density);
                thismod.dz.push_back(dzin);
            }
            din.close();
            this->mods.insert(pair<string,LayeredModel>(sta,thismod));
            thismod.alpha.clear();
            thismod.beta.clear();
            thismod.rho.clear();
            thismod.dz.clear();
        }
        /* manually build the hard coded default model.
        This code tacitly assumes constants like this are correctly 
        cast automatically to double */
        default_model.alpha.push_back(6.0);   default_model.alpha.push_back(8.0);
        default_model.beta.push_back(3.4641);   default_model.beta.push_back(4.6188);
        default_model.rho.push_back(2.69);   default_model.rho.push_back(3.33);
        default_model.dz.push_back(35.0);   default_model.dz.push_back(0.0);
        default_model.name=string("default");
    } catch(...){throw;};
}
StaVariableLayeredSynthetic::StaVariableLayeredSynthetic(string modfile,
        double tsig,double wlev,bool MakeRF,string format)
{
    string base_error("StaVariableLayeredSynthetic single file constructor:  ");
    /* for now only one format allowed */
    if(format!=default_model_format) throw SeisppError(base_error
            + "coding error.  Only accept default format = "
            + default_model_format);
    tsigma=static_cast<float>(tsig);
    wlevel=static_cast<float>(wlev);
    ConvertToRF=MakeRF;
    /* This mixed C and C++ io is less than elegant, but my choice*/
    FILE *fp;
    fp=fopen(modfile.c_str(),"r");
    if(fp==NULL) throw SeisppError(base_error
            + "fopen failure for file="+modfile);
    try {
        LayeredModel thismod;
        char staproto[64];
        int nlayers;
        int i;
        while(fscanf(fp,"%s%d",staproto,&nlayers)!=EOF)
        {
            double Pvin,Svin,rhoin,dzin;
            string sta(staproto);
            for(i=0;i<nlayers;++i)
            {
                fscanf(fp,"%lf%lf%lf%lf",&Pvin,&Svin,&rhoin,&dzin);
                thismod.alpha.push_back(Pvin);
                thismod.beta.push_back(Svin);
                thismod.rho.push_back(rhoin);
                thismod.dz.push_back(dzin);
            }
            this->mods.insert(pair<string,LayeredModel>(sta,thismod));
            thismod.alpha.clear();
            thismod.beta.clear();
            thismod.rho.clear();
            thismod.dz.clear();
        }
        /* WARNING:  these need to be kept consistent with above. Awkward
           potential maintenance issue */
        default_model.alpha.push_back(6.0);   default_model.alpha.push_back(8.0);
        default_model.beta.push_back(3.4641);   default_model.beta.push_back(4.6188);
        default_model.rho.push_back(2.69);   default_model.rho.push_back(3.33);
        default_model.dz.push_back(35.0);   default_model.dz.push_back(0.0);
        default_model.name=string("default");
    }catch(...){throw;};
}
/* Used to find the P arrival time on the P synthetic.   Looks for largest amplitude 
   so this assumes the input is synthetic for z.  Returns the sample number (first 
   sample 0 reference) of peak position. There may be a C++ standard algorithm for
   this, but this is so simple will simply code it here.*/
int find_first_peak(float *d,int n)
{
    int imax=0;
    float dmax=d[0];
    int i;
    for(i=1;i<n;++i)
    {
        if(d[i]>dmax)
        {
            imax=i;
            dmax=d[i];
        }
    }
    return(imax);
}
ThreeComponentSeismogram StaVariableLayeredSynthetic::Compute3C(
        const ThreeComponentSeismogram& parent,Hypocenter& hypo,
        double rlat, double rlon, double relev,string units)
{
    const string base_error("StaVariableLayeredSynthetic::Compute3C method:  ");
    try {
        ThreeComponentSeismogram result(parent);
        result.u.zero();
        /* Velocity models are indexed by sta name so first we have
           to fetch that name.  As usual in seispp the frozen name 
           for the attibute is problematic*/
        string sta=result.get_string("sta");
        LayeredModel vmodel=get_model(sta);
        int nlyrs=vmodel.alpha.size();
//DEBUG
/*
cout << "model for sta="<<sta<<endl;
for(int im=0;im<nlyrs;++im)
{
cout << vmodel.alpha[im] <<" "
  << vmodel.beta[im]<<" "
  << vmodel.rho[im]<<" "
  << vmodel.dz[im]<<endl;
}
*/
        double backaz,azimuth;
        backaz=hypo.seaz(rlat,rlon);
        //dist returns backazimuth, convert to propagation azimuth
        azimuth=backaz+M_PI; 
        /* shift to relative time if necessary */
        double patime;
        if(result.tref == absolute)
        {
            patime=hypo.ptime(rlat,rlon,relev);
            result.ator(patime);
        }
        /* Need slowness vector for incident slowness for the synthetic*/
        SlownessVector uvector(hypo.pslow(rlat,rlon,relev));
        float slow=static_cast<float>(uvector.mag());
        int data_nsamp=result.u.columns();
        /* Kennett's fortran synthetic program uses a power of 2 fft
           so we need to find the next larger power of 2 for size */
        double log2ns=log2(static_cast<double>(data_nsamp));
        // A rather eloborate constuct to remove fraction of log2ns and add one
        int nsamp=static_cast<int>(pow(2.0,static_cast<double>(static_cast<int>(log2ns+1.0))));
//DEBUG
//cout << "nsamp="<<nsamp<<endl;
        float *zp,*rp,*zs,*rs,*trans,*rfr;
	// Intentionally add 4 extra slots at end because of a potential off by one
	// problem created by need for nasty opaque pointers casting float array to a 
	// fortran 32 bit complex
        zp=new float[nsamp+4];
        rp=new float[nsamp+4];
        zs=new float[nsamp+4];
        rs=new float[nsamp+4];
        trans=new float[nsamp+4];
        rfr=new float[nsamp+4];
        int ierr;
        /* Now call the fortran routine.  Because it is fortran
           all arguments have to be passed as pointers.  We use
           the trick that for a vector container we can get a
           array pointer as used here. */
        int fullrsp(1);
        float fdt=static_cast<float>(result.dt);
        kntsyn_(&nlyrs,&(vmodel.alpha[0]),&(vmodel.beta[0]),
                &(vmodel.rho[0]),&(vmodel.dz[0]),&slow,&fullrsp,
                &fdt,&nsamp,&(this->tsigma),&(this->wlevel),
                zp,rp,zs,rs,trans,rfr,&ierr);
        if(ierr)
        {
            delete [] zp;
            delete [] rp;
            delete [] zs;
            delete [] rs;
            delete [] trans;
            delete [] rfr;
            if(ierr==-1) throw SeisppError(base_error
                    + "velocity model for station "
                    + sta + " has more layers than maximum allowed");
            else if(ierr==(-2))
                throw SeisppError(base_error 
                        + "required power of 2 for ns input violated.  Coding error needs to be fixed.");
            else
                throw SeisppError(base_error
                        + " unexpected error code returned by kntsyn");
        }
        double tzerolag(5.0);
        int zerolag;
        if(!ConvertToRF)
        {
            /* this helper hunts for first peak in z component */
            zerolag=find_first_peak(zp,nsamp);
            tzerolag=static_cast<double>(zerolag)*result.dt;
        }
        /* We need to compute the offset in work space from sample 0.
           Since seispp library allows negative t0 this has to be computed
           as here.  Sign oddity because positive lag is positive time shift*/
        int lag0=SEISPP::nint(-(result.t0+tzerolag)/result.dt);
        int ns_to_copy;
        ns_to_copy=result.ns-lag0;
        if(ns_to_copy>nsamp) ns_to_copy=nsamp;
        /* Component 1 is radial and 2 is vertial in internal TRL coordinates.  
           this is reason for indexing in loops below.  The use of a
	   a conditional to avoid negative indices is a bit ineffiicent
           but will leave it unless this proves problematic. */
        int i,ii;
        if(ConvertToRF)
        {
            for(i=0,ii=lag0;i<ns_to_copy;++i,++ii) 
		if(ii>=0) result.u(1,ii)=rfr[i];
        }
        else
        {
            /* This might be faster with the blas, but more confusing
               so use this clearer construct intentionally*/
            for(i=0,ii=lag0;i<ns_to_copy;++i,++ii)
            {
		if(ii>=0)
		{
                    result.u(1,ii)=rp[i];
                    result.u(2,ii)=zp[i];
		}
            }
        }
        /* Done with these work spaces so we can clear them now. */
        delete [] zp;
        delete [] rp;
        delete [] zs;
        delete [] rs;
        delete [] trans;
        delete [] rfr;
        /* We need to construct the transformation matrix since the plane wave
           algorithm computes vertical and radial component data.  We set the
           transformation matrix and call the method of the ThreeComponentSeismogram
           that returns the data to cardinal directions. */
        double a,b;
        a=cos(azimuth);
        b=sin(azimuth);
        result.tmatrix[0][0] = a;
        result.tmatrix[1][0] = b;
        result.tmatrix[2][0] = 0.0;
        result.tmatrix[0][1] = -b;
        result.tmatrix[1][1] = a;
        result.tmatrix[2][1] = 0.0;
        result.tmatrix[0][2] = 0.0;
        result.tmatrix[1][2] = 0.0;
        result.tmatrix[2][2] = 1.0;
        // Perhaps unnecessary, but best to be sure these are set
        result.components_are_cardinal=false;
        result.components_are_orthogonal=true;
        result.rotate_to_standard();
        return(result);
    }catch(...){throw;};
}
LayeredModel StaVariableLayeredSynthetic::get_model(string staname)
{
    map<string,LayeredModel>::iterator it;
    it=mods.find(staname);
    if(it==mods.end())
        return(default_model);
    else
        return((*it).second);
}
TimeSeries StaVariableLayeredSynthetic::ComputeScalar(const TimeSeries& parent,
        Hypocenter& hypo,double rlat, double rlon, double relev,string type)
{
    cout << "StaVariableLayeredSynthetic::ComputeScalar procedure not implemented"<<endl;
    exit(-1);
}
TimeSeries StaVariableLayeredSynthetic::ComputeScalar(int nsamp,
        double dt, Hypocenter& hypo,
        double rlat, double rlon, double relev,string type)
{
    // This method can build a parent TimeSeries object, call the above, and
    // return it.  
    cout << "StaVariableLayeredSynthetic::ComputeScalar plain procedure not implemented"<<endl;
    exit(-1);
}
ThreeComponentSeismogram StaVariableLayeredSynthetic::Compute3C(
        int nsamp, double dt,
        Hypocenter& hypo, double rlat, double rlon, double relev,
        string options)
{
    try {
        ThreeComponentSeismogram pattern(nsamp);
        pattern.live=true;
        stringstream ss(options);
        double t0;
        string sta;
        /* Far from elegant pushing this in a string with a required order */
        ss >> sta;
        ss >> t0;
        /* if the t0 string is missing the above sets it 0 so we can default this way. */
        if(t0==0.0) t0=-10.0;
        pattern.dt=dt;
        pattern.t0=t0;
        pattern.put("sta",sta);
        return(this->Compute3C(pattern,
                    hypo,rlat,rlon,relev,options));
    }catch(...){throw;};
}
/* Could put this inline but done this way as I'm not positive this simple code will work. */
void StaVariableLayeredSynthetic::set_model(string sta, LayeredModel& m)
{
    mods[sta]=m;
}
