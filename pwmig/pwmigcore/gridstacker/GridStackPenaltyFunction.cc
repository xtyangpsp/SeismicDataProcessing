#include "perf.h"
#include "Metadata.h"
#include "GridStackPenaltyFunction.h"
GridStackPenaltyFunction::GridStackPenaltyFunction(Metadata& md)
{
   try {
	floor=md.get_double("robust_weight_floor");
	string pfuncname=md.get_string("PenaltyFunction");
	if(pfuncname=="simple_coherence")
	{
		gpt=COH;
		cohpow=1.0;
	}
	else if(pfuncname=="coherence_power_function")
	{
		gpt=COHPOW;
		cohpow=md.get_double("power_factor");
	}
	else if(pfuncname=="dbxcor_loss_function")
	{
		gpt=DBXCOR;
		cohpow=md.get_double("power_factor");
	}
	else
	{
		cerr << "GridStackPenaltyFunction constructor:  "
			<< "Unknown name for PenaltyFuntion="
			<< pfuncname<<endl;
		exit(-1);
	}
	
    }
    catch (MetadataGetError mderr)
    {
	mderr.log_error();
	exit(-1);
    }
}
double GridStackPenaltyFunction::weight(int nd, double *d, GCLvectorfield3d& d0)
{
	double sumsqr,sumsqd,r;
	int i,j,k,l;
	int ii;
	double coh,wt;
	int ntest=(d0.n1)*(d0.n2)*(d0.n3)*(d0.nv);
	double ddotd0;  // needed in DBXCOR version
	/* This should be converted to an exception if this code is ever exported outside
	this program.  Here this is just a sanity check for coding errors.  Faster and easier
	to just trap this with an exit. */
	if(ntest!=nd)
	{
		cerr << "GridStackPenaltyFunction::weight method:  coding error"<<endl;
		cerr << "nd="<<nd<<" but product grid dimensions="<<ntest<<endl;
		exit(-1);
	}
	/* WARNING:  this is heavily dependent on the internal structure of the field 
	array used in the current GCLgrid library. It assumes a precise index order
	mapping to a contiguous vector of doubles.  */
	ii=0;
	sumsqr=0.0;  sumsqd=0.0;  ddotd0=0.0;  
	for(i=0;i<d0.n1;++i)
	  for(j=0;j<d0.n2;++j)
	    for(k=0;k<d0.n3;++k,ii+=5)
	    {
		/* this is assumed to have been set for d when it was created from a grid file */
		if( (d[ii+4]>0.0) && (d0.val[i][j][k][4]>0.1) )
		{
//DEBUG
//cout << i <<" "<<j<<" "<<k;
//cout << "   d[ii+4]="<<d[ii+4]<<" d0.val4="<<d0.val[i][j][k][4]<<endl;
			for(l=0;l<3;++l)
			{
//DEBUG
/*
int itest;
itest=i*d0.n2*d0.n3*d0.nv + j*d0.n3*d0.nv + k*d0.nv + l;
cout << "itest="<<itest<<" ii+l="<<ii+l<<endl;
*/
				r=d[ii+l]-d0.val[i][j][k][l];
//cout <<"r="<<r<< " = d[ii+1]="<<d[ii+l] <<" - d0val="<<d0.val[i][j][k][l]<<endl;
				sumsqr+=r*r;
				sumsqd+=d[ii+l]*d[ii+l];
				/* This is needed only in dbxcor weight scheme, but easier to cost to 
				compute it here is tiny compared to overhead (an potential error) 
				in writing a separate loop */
				ddotd0+=d[ii+l]*d0.val[i][j][k][l];
			}
		}
	    }

	if(sumsqr>=sumsqd)
		coh=0.0;
	else
		coh=1.0-sqrt(sumsqr/sumsqd);
	switch(gpt)
	{
	case COH:
		wt=coh;
		break;
	case DBXCOR:
		wt=ddotd0/sqrt(sumsqd*sumsqr);
		wt=pow(wt,cohpow);
		break;
	case COHPOW:
	default:
		wt=pow(coh,cohpow);
		break;
	}
	if(wt<floor) wt=floor;
	return(wt);
}
