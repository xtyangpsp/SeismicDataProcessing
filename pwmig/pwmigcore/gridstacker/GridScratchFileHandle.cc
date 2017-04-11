#include <stdio.h>
#include <list>
#include "seispp.h"
#include "GridScratchFileHandle.h"
using namespace std;
using namespace SEISPP;

MemberGrid::MemberGrid(const MemberGrid& parent)
{
	gridname=parent.gridname;
	fieldname=parent.fieldname;
	baseweight=parent.baseweight;
	reswt=parent.reswt;
}

MemberGrid& MemberGrid::operator= (const MemberGrid& parent)
{
	if(this!=&parent)
	{
		gridname=parent.gridname;
		fieldname=parent.fieldname;
		baseweight=parent.baseweight;
		reswt=parent.reswt;
	}
	return(*this);
}

/* More rigorous test for equality than operator== required for this code */
bool compare_equal(GCLvectorfield3d& g1, GCLvectorfield3d& g2)
{
	if(g1!=g2) return(false);
	if(g1.n1!=g2.n1) return(false);
	if(g1.n2!=g2.n2) return(false);
	if(g1.n3!=g2.n3) return(false);
	if(g1.nv!=g2.nv) return(false);
	return(true);
}
const double null_field_value(-1.0);  // inserted in solid angle field (3) when data are considered null

GridScratchFileHandle::GridScratchFileHandle(GCLvectorfield3d& mastergrid,
		list<MemberGrid>& mgl, bool normalize, double cutoff)
{
	const string base_error("GridScratchFileHandle constructor:  ");
	fp=tmpfile();
	if(fp==NULL) throw SeisppError(base_error
			+ string("tmpfile() failed"));
	n1=mastergrid.n1;
	n2=mastergrid.n2;
	n3=mastergrid.n3;
	nv=mastergrid.nv;
	nval=n1*n2*n3*nv;
	current_member=0;
	nmembers=0;  // Could use size on list, but better to count grids loaded
	list<MemberGrid>::iterator mptr;
	for(mptr=mgl.begin();mptr!=mgl.end();++mptr)
	{
	    try{
                GCLvectorfield3d current(mptr->fieldname);
		if(normalize) 
		{
			// This should perhaps be a function, but will inline it for efficiency 
			int i,j,k,l;
			for(i=0;i<current.n1;++i)
			 for(j=0;j<current.n2;++j)
			  for(k=0;k<current.n3;++k)
			  {
				if(current.val[i][j][k][3]>cutoff)
				{
					for(l=0;l<3;++l) current.val[i][j][k][l]/=current.val[i][j][k][3];
				}
				else
				{
					current.val[i][j][k][4]=-1.0;
					/* We zero data below cutoff to allow a vector stack of data.  Then
					we deal with normalization by counting points with field 3 positive.  
					see stacking code */
					for(l=0;l<3;++l) current.val[i][j][k][l]=0.0;
				}
			   }
		}
		if(compare_equal(mastergrid,current))
		{
//DEBUG
vector<double> test;
int i,j,k,l;
for(i=0;i<current.n1;++i)
for(j=0;j<current.n2;++j)
for(k=0;k<current.n3;++k)
test.push_back(current.val[i][j][k][2]);
			if(current.nv==mastergrid.nv)
				fwrite(current.val[0][0][0],sizeof(double),nval,fp);
			else
			{
				cerr << base_error
					<< "vector field size mismatch for fieldname="
					<< mptr->fieldname
					<< endl
					<< "Required nv of master="
					<<mastergrid.nv
					<< "this member's nv="
					<< current.nv
					<< endl
					<< "Skipping data for this grid"
					<< endl;
			}
		}
		else
		{
			GCLvectorfield3d scratch(mastergrid);
			scratch.zero();
			scratch += current;
			fwrite(scratch.val[0][0][0],sizeof(double),nval,fp);
		}
		++nmembers;
		++current_member;
	    } catch (int ierr) {
		cerr << base_error
			<< "GCLvectorfield3d constructor failed for grid="
			<<mptr->gridname
			<< " fieldname="
			<< mptr->fieldname
			<<endl
			<<"These data will be dropped from the stack"
			<<endl;
	    }    
	}
}
		
GridScratchFileHandle::~GridScratchFileHandle()
{
	if(fp!=NULL) fclose(fp);
}

void GridScratchFileHandle::rewind()
{
	int iret=fseek(fp,0L,SEEK_SET);
	current_member=0;
}
/*
// The format of this is unquestionably wrong.  NO references to consult
void GridScratchFileHandle::operator ++()
{
	double *buffer=new double[nval];
	int nread=fread(buffer,sizeof(double),nval,fp);
	delete [] buffer;
	if(nread!=nval) throw SeisppError(
		string("GridScratchFileHandle::operator++:  fread error on scratch file"));
	++current_member;

}
*/
int GridScratchFileHandle::load_next(double *buffer)
{
	int nread=fread(buffer,sizeof(double),nval,fp);
	if(nread!=nval) throw SeisppError(
		string("GridScratchFileHandle::load_next:  fread error"));
	++current_member;
	return(nread);
}
bool GridScratchFileHandle::at_eof()
{
	if(current_member==nmembers)
		return true;
	else
		return false;
}

