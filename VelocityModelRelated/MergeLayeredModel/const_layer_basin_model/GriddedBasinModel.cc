/* Takes a list of isopach surfaces constructed with isopach2gcl and builds a
   constant velocity and density basin model from these.   Vp, Vs, and density
   of each unit are provided as input */
#include <string>
#include <vector>
#include <sstream>
#include "gclgrid.h"
#include "GriddedBasinModel.h"
#include "SeisppError.h"
using namespace std;
using namespace SEISPP;
/* This assumes the first member grid geometry is same as others.  It also
assumes the order starts at the surface and is ordered downward AND the 
tabulated z values are elevation of the surface. */
GriddedBasinModel::GriddedBasinModel(vector<GCLscalarfield> s,
        vector<double> vp, vector<double> vs, vector<double>rho,
        string nm)
        : GCLvectorfield3d(s[0].n1,s[0].n2,s.size(),4)
{
    string base_error("GriddedBasinModel constructor:  ");
    int nsurfaces=s.size();
    if(nsurfaces!=vp.size())
    {
        throw SeisppError(base_error
            + "size mismatch between list of isopach files and vp values");
    }
    if(nsurfaces!=vs.size())
    {
        throw SeisppError(base_error
            + "size mismatch between list of isopach files and vs values");
    }
    if(nsurfaces!=rho.size())
    {
        throw SeisppError(base_error
            + "size mismatch between list of isopach files and density values");
    }
    int i,j,k,kk,l;
    this->name=nm;
    model_name=nm;
    /* This is ugly as I couldn't find a way to do this through the
       constructor line.  This is copying all public attributes defined
       in BasicGCLgrid.  This is a flaw in the gclgrid library design */
    lat0=s[0].lat0;
    lon0=s[0].lon0;
    r0=s[0].r0;
    azimuth_y=s[0].azimuth_y;
    dx1_nom=s[0].dx1_nom;
    dx2_nom=s[0].lat0;
    dx3_nom=0.2;  // arbitrary but ok for Ill basin
    /* These will be wrong but I don't think it matters */
    i0=s[0].i0;
    j0=s[0].j0;
    k0=0;
    for(i=0;i<3;++i) translation_vector[i]=s[0].translation_vector[i];
    for(i=0;i<3;++i)
        for(j=0;j<3;++j) gtoc_rmatrix[i][j]=s[0].gtoc_rmatrix[i][j];
    /* end of stuff that really belongs in parent constructor.
       First transfer the coordinates of the top surface.   These
       data are elevations so this datum is not sea level */
    int nz=s.size();
    
    for(k=0,kk=nz-1;kk>=0;++k,--kk)
    {
        /* GCLgrids require the cells be right handed for lookup member
           to work correctly.  Hence, we have to reverse the k index in
           the grid compared to the input.  This pair of backward running
        indices is intrinsically dangerous. */
        // Sanity check on sizes is the only check we make
        // Not as robust as it could be but this is not intended for
        // wide release
        if((s[k].n1 != this->n1) || (s[k].n2 != this->n2) )
            throw SeisppError(base_error + "size mismatch of n1 and n2 in member");
        /* Surfaces for each level are already in the right 
           coordinates so we can just copy them. */
        for(i=0;i<n1;++i)
            for(j=0;j<n2;++j)
            {
                x1[i][j][kk]=s[k].x1[i][j];
                x2[i][j][kk]=s[k].x2[i][j];
                x3[i][j][kk]=s[k].x3[i][j];
                
                this->val[i][j][kk][0]=vp[k];
                this->val[i][j][kk][1]=vs[k];
                this->val[i][j][kk][2]=rho[k];
                
                /* zero the thickness field for now so we can 
                   run the sanity check that follow */
                this->val[i][j][kk][3]=0.0;
                
                //this->val[i][j][kk][3]=s[k].val[i][j];
                //cout<<"thickness= "<<s[k].val[i][j]<<endl;
            }
    }
    /* Now fill in the values for the vector field.   In all cases we 
       write the field values in vector slots with 0-vp, 1-vs, 2-rho,
       and 3-thickness.  Eventually these could be space variable. */
    /*
    for(i=0;i<n1;++i)
        for(j=0;j<n2;++j)
            for(k=0,kk=n3-1;k<n3;++k,--kk)
            {
                this->val[i][j][kk][0]=vp[k];
                this->val[i][j][kk][1]=vs[k];
                this->val[i][j][kk][2]=rho[k];
                */
                /* zero the thickness field for now so we can 
                   run the sanity check that follow */
                //this->val[i][j][kk][3]=0.0;
                /*
                this->val[i][j][kk][3]=s[k].val[i][j];
                cout<<"thickness= "<<s[k].val[i][j]<<endl;
            }
    */
    /* The grid geometry used here will fail if any layer has zero thickness.
       To reduce that problem this section forces a minimum thickness of 
       any horizon to the const, frozen value.  Note the actual files
       have a 2 m minimum, so this is somewhat redundant.*/
    const double MinThickness(0.001);  // This is 1 m 
    for(i=0;i<n1;++i)
        for(j=0;j<n2;++j)
        {
            double zk(0.0),zkm,currentthickness;
            int sanitycounter;
            const int insane(1000);
            //debugging
    		//vector<double> z;
    		//z.push_back(0.0); //surface depth or top depth of the first layer.
            for(k=n3-1;k>0;--k)
            //for(k=n3-2;k>=0;--k)
            {
                //cout<<"Working on layer = "<<k<<endl;
                /*
                currentthickness=this->val[i][j][k][3];
                zk += currentthickness;
                //z.push_back(zk);
                
                Geographic_point gp=this->geo_coordinates(i,j,k);
				gp.r=r0_ellipse(gp.lat) - zk;
				Cartesian_point cp=this->gtoc(gp);
				this->x1[i][j][k]=cp.x1;
				this->x2[i][j][k]=cp.x2;
				this->x3[i][j][k]=cp.x3;
                */
                
                /* remember order in k index starts at base and works to top */
                zk=this->depth(i,j,k);
                zkm=this->depth(i,j,k-1); //zkm is zk minus
                //cout<<"zk= "<<zk<<", zkm= "<<zkm<<endl;
                /*
                zk=this->depth(i,j,k+1);
                zkm=this->depth(i,j,k); //zkm is zk minus
                */
                
                if(zk>=zkm)
                {
                    double zk0=zk;
                    zk+=MinThickness;
                    sanitycounter=0;
                    while((zk>=zkm) && (sanitycounter<insane))
                    {
                        //debug
                        //cerr << "zk="<<zk<<endl;
                        zk -= MinThickness;
                        ++sanitycounter;
                    }
                    //debug
                    //cout<<"insane="<<insane<<", sanitycounter = "<<sanitycounter<<endl;
                    if(sanitycounter>=insane) throw SeisppError(base_error
                            + "Absurd data input.  Layers likely out of order.");
                    cerr << "Adjusted depth of horizon "<<k<<" from "
                        << zk0 <<" to " << zk<<endl;
                    Geographic_point gp=this->geo_coordinates(i,j,k);
                    gp.r=r0_ellipse(gp.lat) - zk;
                    Cartesian_point cp=this->gtoc(gp);
                    this->x1[i][j][k]=cp.x1;
                    this->x2[i][j][k]=cp.x2;
                    this->x3[i][j][k]=cp.x3;
                }
                
                //this->val[i][j][k+1][3]=zkm-zk; // odd indexing, but an evolutionary result
                this->val[i][j][k][3]=zkm-zk; 
            }
        }

    //This may not be necessary but a simple way to set bounding attributes
    this->compute_extents();
}
GriddedBasinModel::GriddedBasinModel(string fname) 
    : GCLvectorfield3d(fname)
{
    /* With default format the file name base is the field name */
    model_name=fname;
}
void GriddedBasinModel::save(string fname, string dir,string format)
{
    /* We downcast this to the parent and save it that way.   This allows
       the constructor to be little more than an alias for a GCLvectorfield3d.
       It is, however, somewhat bad form and not something to build on.
       Doing this because this is not intended to be a transportable code
       but one that is more or less a prototype for a more general version
       I would ultimately release. */
    /* that does not work - typeid.name seems to always return the most
       derived class. Hence the real real bubba solution following.
    GCLvectorfield3d *parent;
    parent=dynamic_cast<GCLvectorfield3d *>(this);
    parent->save(fname,dir,format);
    */
    GCLvectorfield3d parent;
    parent=dynamic_cast<GCLvectorfield3d &>(*this);
    parent.save(fname,dir,format);
}
LayeredModel GriddedBasinModel::fetch_model(double lat0, double lon0)
{
    try {
        const string base_error("GriddedBasinModel::fetch_model:  ");
        int iret=this->lookup(lat0,lon0,r0_ellipse(lat0)-1.0);
        if(iret<0) 
        {
            stringstream ss;
            ss << base_error << "lookup failed for point at latitude "
                << deg(lat0) <<" and longitude "<<deg(lon0)<<endl;
            throw SeisppError(ss.str());
        }
        /* for now we simply used the index of the lookup point.  More
           than adequate for use in OIINK  and Illinois Basin model */
        int pos[3];
        this->get_index(pos);
        LayeredModel result;
        stringstream ss;
        ss << this->name<<"_"<<pos[0]<<"_"<<pos[1];
        result.name=ss.str();
        for(int k=0;k<n3;++k)
        {
            double alpha,beta,rho,dz;
            // n3-1 because we reference top surface always 
            alpha=static_cast<float>(this->val[pos[0]][pos[1]][n3-1][0]);
            beta=static_cast<float>(this->val[pos[0]][pos[1]][n3-1][1]);
            rho=static_cast<float>(this->val[pos[0]][pos[1]][n3-1][2]);
            dz=static_cast<float>(this->val[pos[0]][pos[1]][n3-1][3]);
            result.alpha.push_back(alpha);
            result.beta.push_back(beta);
            result.rho.push_back(rho);
            result.dz.push_back(dz);
        }
        return(result);
    }catch(...){throw;};
}
        
