#include <float.h>
#include "coords.h"
#include "perf.h"
#include "dmatrix.h"
#include "gclgrid.h"
#include "seispp.h"
using namespace std;
using namespace SEISPP;
#include "pwmig.h"

// Helper function needed below

/* computes the set of local
verticals along the ray path in the GCLgrid cartesian system.  Because local
vertical is a geographic concept the parent GCLGrid, g, needs to be passed
to get geo coordinates.    The g grid is assumed to be consistent with 
the coordinates of path.
*/
dmatrix *compute_local_verticals(GCLgrid& g, dmatrix& path)
{
	int *sz,nx;
	Geographic_point x0_geo;
	Cartesian_point x0_c;
	const double DR=100.0;

	double dx[3];
	double nrmdx;
	int i,j;

	sz = path.size();
	// assume number rows = 3.   could test and throw exception, but this
	// function is not viewed as a library function
	nx=sz[1];
	delete [] sz;
	dmatrix *vptr = new dmatrix(3,nx);
	dmatrix& verticals=*vptr;
	for(i=0;i<nx;++i)
	{
		x0_geo = g.ctog(path(0,0),path(1,0),path(2,0));
		x0_geo.r -= DR;
		x0_c = g.gtoc(x0_geo);
		verticals(0,i) = (path(0,i) - x0_c.x1)/DR;
		verticals(1,i) = (path(1,i) - x0_c.x2)/DR;
		verticals(2,i) = (path(2,i) - x0_c.x3)/DR;
	}
	return(vptr);
}

/* computes depth dependent transformation matrix assuming a scattering
is specular.  That is, transformation matrices will yield a form of
T,R,L as output.
*/


dmatrix Ray_Transformation_Operator::apply(dmatrix& in)
{
	int i,j;
	double *p;
	int nrow=in.rows();
	int ncol=in.columns();
	// This could be trapped as an exception but since this
	// is expected to only be used in this program this is
	// simpler
	//Old form:  if( (nrow!=3) || (ncol!=npoints) )
	if( (nrow!=3) || (ncol>npoints))
	{
		cerr << "Coding error:  Ray_Transformation_Operator has "
			<< npoints << "elements but matrix passed is "
			<< nrow << "X" << ncol << endl
			<< "Data matrix samples (columns) must be smaller "
			<< "than number of operator elements"<<endl;
		exit(-1);
	}
	dmatrix out(nrow,ncol);

	// 
	// This probably could be replaced with matrix vector
	// multiples, but I took this approach under the (perhaps
	// invalid) assumption this would be faster since it 
	// uses pointer arithmetic avoiding all the operator()
	// calls that would otherwise be needed.
	//
	for(i=0;i<ncol;++i)
	{
		double work[3];
		work[0]=in(0,i);
		work[1]=in(1,i);
		work[2]=in(2,i);
		// This is a multiply of U*d
		p=U[i].get_address(0,0);
		out(0,i)=p[0]*work[0]+p[3]*work[1]+p[6]*work[2];
		out(1,i)=p[1]*work[0]+p[4]*work[1]+p[7]*work[2];
		out(2,i)=p[2]*work[0]+p[5]*work[1]+p[8]*work[2];
	}
	return(out);
}
// Simple constuctor 
Ray_Transformation_Operator::Ray_Transformation_Operator(int np)
{
	npoints = np;
        for(int i=0;i<np;++i)
        {
                dmatrix *dmtp = new dmatrix(3,3);
                U.push_back(*dmtp);
		delete dmtp;
        }
}
/* Constructor for simple case with all matrices going to surface R,T,L coordinates.
The primary contents is a vector of dmatrix object holding (in this case) the
same transformation matrix at every point.  

Arguments:
	g - reference GCLgrid object.  Mainly used here as the coordinate reference
		frame for path
	path - 3xNp matrix of points that defines the ray path
	azimuth - angle (radians) of the ray propagation direction 
		at the surface.  Note this is geographical azimuth NOT the
		phi angle in spherical coordinates.
*/
Ray_Transformation_Operator::Ray_Transformation_Operator(GCLgrid& g, 
	dmatrix& path,double azimuth)
{
	int np=path.columns();
	dmatrix U0(3,3);  // final matrix copied to all elements here
	double x[3];  //work vector to compute unit vectors loaded to U0
	Geographic_point x0_geo;
	Cartesian_point x0_c;
	double x_vertical[3];  // point to local vertical direction at x_geo
	const double DR=100.0;
	double xdotxv,theta;
	double a,b,c,d;
	int i;

	//first we need the geographic coordinates of the ray emergence point
	x0_geo = g.ctog(path(0,0),path(1,0),path(2,0));
	x0_geo.r -= DR;
	x0_c = g.gtoc(x0_geo);
	x_vertical[0] = (path(0,0) - x0_c.x1)/DR;
	x_vertical[1] = (path(1,0) - x0_c.x2)/DR;
	x_vertical[2] = (path(2,0) - x0_c.x3)/DR;
	// Get the L direction.  Loop is used to make sure the
	// distance is not tiny which could lead to an NaN in
	// the normalized vector.  Size assumes we are using km
	// as the distance scale.  
	const double dxmin(2.0);
	double nrmx;
	for(i=0;i<path.columns();++i)
	{
		x[0]= path(0,0) - path(0,i);
		x[1]= path(1,0) - path(1,i);
		x[2]= path(2,0) - path(2,i);
		nrmx=dnrm2(3,x,1);
		if(nrmx>dxmin) break;
	}
	dscal(3,1.0/nrmx,x,1);  // normalize to unit vector
	// need to check for a vertical vector and revert to identity 
	// to avoid roundoff problems
	xdotxv = ddot(3,x,1,x_vertical,1);
	/* conservative test for vertical incidence*/
	a=cos(azimuth);
	b=sin(azimuth);
	if((1.0-fabs(xdotxv))<=FLT_EPSILON)
	{
		// L vertical means theta=0.0
		c=1.0;
		d=0.0;
	}
	else
	{
		// dot product is preserved so we can use this here
		theta = acos(xdotxv);
	        c = cos(theta);
	        d = sin(theta);
	}
	/* This was wrong .  Retained temporarily to understand some
	incorrect older results.
	U0(0,0) = a*c;
        U0(1,0) = -b;
        U0(2,0) = a*d;
        U0(0,1) = b*c;
        U0(1,1) = a;
        U0(2,1) = b*d;
        U0(0,2) = -d;
        U0(1,2) = 0.0;
        U0(2,2) = c;
	*/
	/* new version */
	U0(0,0) = a;
        U0(1,0) = b*c;
        U0(2,0) = b*d;
        U0(0,1) = -b;
        U0(1,1) = a*c;
        U0(2,1) = a*d;
        U0(0,2) = 0.0;
        U0(1,2) = -d;
	U0(2,2) = c;
	
	for(int i=0;i<np;++i)
	{
		// This works because container push_back calls copy constructor on *dtmp
		dmatrix *dtmp = new dmatrix(3,3);
		*dtmp = U0;
		U.push_back(*dtmp);  
		delete dtmp;
	}
	npoints=U.size();
}
/* Now the much more complex algorithm for the case where we assume the ray path
is a specular reflection.  Then we have a different ray transformation matrix 
at each point.  We compute this using this algorithm:
  The L coordinate turns upward.  At the scatter point we get it from the scattered ray path
  direction and assume that is one thing we can anchor on.

  The scattered S has a polarization in the plane of scattering formed by the incident P and
  scattered S.  The complementary component is then readily found by the cross product of the
  P and S paths.  

  We get the scattered S direction (at the scattering point) from a cross product between 
  local L and the local "SH-like" direction derived in the previous step.

  We propagate the components to the surface with a simple directional change.  L remains
  L  and the S components obey Snell's law rules.  

Arguments:
	g-parent GCLgrid that defines coordinate system (only used for coordinate
		transformations here)
	path - 3xnp matrix of points defining the path.  It is assumed path starts at
		the surface and is oriented downward.
	azimuth - angle (radians) of the ray propagation direction 
		at the surface.  Note this is geographical azimuth NOT the
		phi angle in spherical coordinates.
	gamma_P - 3xnp matrix of tangent (unit) vectors for incident wave ray path
		at each point in path.  It is assumed gamma_P vectors point along
		the direction of propagation of the incident wavefield (nominally upward)

The output array of matrices when applied to data will yield output with x1 = generalized
T, x2= generalized R, and x3=L.  

*/
Ray_Transformation_Operator::Ray_Transformation_Operator(GCLgrid& g, 
	dmatrix& path,double azimuth, dmatrix& gamma_P)
{
	int np=path.columns();
	Geographic_point x0_geo;
	Cartesian_point x0_c;
	double x_vertical[3];  // point to local vertical direction at x_geo
	const double DR=100.0;
	const double PARALLEL_TEST=FLT_EPSILON;  // a conservative test for zero length double vector
	int i,j;

	// First we need to compute a set of ray path tangents
	dmatrix *tanptr;
	tanptr = ray_path_tangent(path);
	dmatrix &tangents = *tanptr;

	// We next need a set of local vertical vectors. 
	dmatrix *vptr;
	vptr = compute_local_verticals(g,path);
	dmatrix &local_verticals = *vptr;

	// We call the simpler constructor first above and use it to provide
	// transformation matrices to ray coordinates -- the starting point here
	// For old C programmers like me:  This will be deleted when it goes out
	// of scope when the function exits so we don't need a delete below.
	Ray_Transformation_Operator Raytrans0(g,path,azimuth);

	for(i=0;i<np;++i)
	{
		double Lscatter[3],Tscatter[3],Rscatter[3];
		double nu0[3];  // unit vector in direction gamma_P
                /* Tp, Rp, and Zp are a nonorthogonal set of vectors used
                   to construct the Lsc, Rsc, and Tsc orthogonal basis below.
                   Zp is local vertical, Tp is tangential in the conventional
                   mode, and Rp is a radial.  Tp will be perpendicular to Zp
                   and Rp, BUT Rp is not perpendicular to Zp but nu0 instead.
                   Kind of a confusing set of symbols, but inherited from
                   some earlier mixups and I don't want to botch the changes.*/
		double Zp[3],Rp[3],Tp[3];  // radial and tangential for S ray path 
		dmatrix work(3,3);
		// copy the tangent vector Lscatter from the tangent vectors 
		// matrix for convenience and clarity.  Assume Lscatter is unit vector
		dcopy(3,tangents.get_address(0,i),1,Lscatter,1);
		// do the same for the incident ray direction but here we want to
		// normalize to a unit vector
		dcopy(3,gamma_P.get_address(0,i),1,nu0,1);
		dr3norm(nu0);
		// First derive the S ray path propagation coordinate system
		// first copy Zp
		dcopy(3,local_verticals.get_address(0,i),1,Zp,1);
		// Now get radial from Lscatter = S ray path tangent.  Have to 
		// handle case with Lscatter vertical carefully.  In that case, 
		// derive radial from nu0.
		dr3cros(Zp,Lscatter,Tp);
		if(dnrm2(3,Tp,1)<PARALLEL_TEST)
                {
			dr3cros(Zp,nu0,Tp);
			// excessively paranoid, but could happen
			if(dnrm2(3,Tp,1)<PARALLEL_TEST)
                        {
                                cerr << "Warning:  cannot handle singular geometry"
					<< endl
				<< "Antipodal event scattered to vertical incidence"
					<< endl
				<< "Using T=x1 and R=x2 of GCL coordinate system"
					<< endl
				<< "Probable discontinuity at normal incidence"<<endl;
				Tp[0]=1.0;  Tp[1]=0.0;  Tp[2]=0.0;
			}
		}
		dr3norm(Tp);
		// Radial is now derived as Lscatter X Tp to make right handed cartesian
		dr3cros(Lscatter,Tp,Rp);
		// Normalization dr3norm(Rp) not necessary here because Tp and Lscatter
		// are orthogonal unit vectors
		//
		// Now we need to derive the specular reflection SV vector (Rscatter)
		// and SH vector (Tscatter).  A confusion is this is much like the
		// way we computed Rp and Tp above but using different vectors as
		// building blocks.  First, Tscatter is derived from cross product of
		// nu0 and Lscatter.  Sign is important.  We aim for a coordinate system
		// where T,R,L = x1,x2,x3 form a right handed coordinate system
		// This happens to be a nonstandard convention (textbooks would
		// use L,R,T = x,y,z) but I'm deriving too much code from 
		// multiwavelet code that used this convention.  
		//
                /* Tscatter is the x1 coordinate.  We do the cross product in
                   this order to be consistent with Tp, Rp, and nu0 that is the
                   set of basis vectors rotating to the T,R,L directions in
                   the conventional way.   An earlier version had this 
                   incorrectly reversed which produces sign issues. */
                dr3cros(nu0,Lscatter,Tscatter);
		// as above, need to handle case when L and T are parallel
		if(dnrm2(3,Tscatter,1)<PARALLEL_TEST)
		{
			// In this case, Tscatter will be the same as Tp
			dcopy(3,Tp,1,Tscatter,1);
		}
		// Do need to normalize
		dr3norm(Tscatter);  
		// In a right handed coordinate system x2 = x3 X x1 = L X T
		dr3cros(Lscatter,Tscatter,Rscatter);

		//
		/* Now we derive the transformation matrix.  The matrixs computed
		computed here will transform a vector in ray coordinates to a vector
		with in the specular reflection coordinate system used here.
		Note that for the present case of a 1d medium this is a simple
		rotation around x3, but we derive it from dot products.  
		In the symbols of this algorithm it transforms from the 
		Tp,Rp, Lscatter basis to that of Tscatter, Rscatter, Lscatter.
		This works because we assume a simple ray theory propagator to
		project the wavefield to the depth of the scattering point.
		(Note: my working theory notes derived the inverse of this
		transform = transpose of this one.
		*/

		work(0,0) = ddot(3,Tscatter,1,Tp,1);
		work(0,1) = ddot(3,Tscatter,1,Rp,1);
		work(0,2) = 0.0;
		work(1,0) = ddot(3,Rscatter,1,Tp,1);
		work(1,1) = ddot(3,Rscatter,1,Rp,1);
		work(1,2) = 0.0;
		work(2,0) = 0.0;
		work(2,1) = 0.0;
		work(2,2) = 1.0;
		U.push_back(work*Raytrans0.U[i]);
	}
	npoints=U.size();

	delete tanptr;
	delete vptr;
}
// Always need a copy constructor.  Not currently used in pwmig, 
// but have learned the hard way a copy constructor is sometimes
// implied.
Ray_Transformation_Operator::Ray_Transformation_Operator(const Ray_Transformation_Operator& pat)
{
	npoints=pat.npoints;
	U=pat.U;
}
