#include "SlownessVectorMatrix.h"
namespace PWMIG {
using namespace SEISPP;
using namespace std;

SlownessVectorMatrix::SlownessVectorMatrix()
{
    nrow=0;
    ncol=0;
}
SlownessVectorMatrix::SlownessVectorMatrix(int nr, int nc)
{
    nrow=nr;
    ncol=nc;
    int n=nr*nc;
    uarray.reserve(n);
    int i;
    for(i=0;i<n;++i) 
        uarray.push_back(SlownessVector(0.0,0.0,0.0));
}
SlownessVectorMatrix::SlownessVectorMatrix(const SlownessVectorMatrix& parent)
                            : uarray(parent.uarray)
{
    nrow=parent.nrow;
    ncol=parent.ncol;
}
SlownessVector& SlownessVectorMatrix::operator()(int i, int j)
{
    if( (i<0) || (i>=nrow) || (j<0) || (j>=ncol) )
    {
        stringstream ss;
        ss << "SlownessVectorMatrix::operator(): row="<<i
            <<" column="<<j<<" is outside bounds of matrix size="
            << nrow <<"X"<<ncol<<endl;
        throw SeisppError(ss.str());
    }
    int n;
    return(this->uarray[i+nrow*j]);
}
SlownessVectorMatrix& SlownessVectorMatrix::operator=
                        (const SlownessVectorMatrix& parent)
{
    if(this!=&parent)
    {
        nrow=parent.nrow;
        ncol=parent.ncol;
        uarray=parent.uarray;
    }
    return *this;
}
/* Made a friend instead of part of object*/
ostream& operator<<(ostream& ostr, SlownessVectorMatrix& svm)
{
    ostr << "row index, column index, ux, uy"<<endl;
    int i,j;
    for(i=0;i<svm.nrow;++i)
        for(j=0;j<svm.ncol;++j)
        {
            SlownessVector uij=svm(i,j);
            ostr << i <<" "
                << j << " "
                << uij.ux << " "
                << uij.uy << endl;
        }
    return ostr;
}
}
