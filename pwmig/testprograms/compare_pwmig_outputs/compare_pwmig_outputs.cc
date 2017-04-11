#include <string>
#include "gclgrid.h"
using namespace std;
void usage()
{
   cerr << "compare_pwmig_files f1 f2 [-o diff]"<<endl
        << " -o writes f1-f2 to diff"<<endl;
}
double L2norm(GCLvectorfield3d& f)
{
    int i,j,k,l;
    double nrm(0.0);
    for(i=0;i<f.n1;++i)
        for(j=0;j<f.n2;++j)
            for(k=0;k<f.n3;++k)
                for(l=0;l<3;++l)
                    nrm+=(f.val[i][j][k][l])*(f.val[i][j][k][l]);
    return sqrt(nrm);
}
int main(int argc, char **argv)
{
    if(argc<3) usage();
    string f1(argv[1]);
    string f2(argv[2]);
    string ofile;
    bool write_output_file(false);
    int i,j,k,l;
    for(i=3;i<argc;++i)
    {
        string sarg(argv[i]);
        if(sarg=="-o")
        {
            ++i;
            if(i>=argc) usage();
            ofile=string(argv[i]);
            write_output_file=true;
        }
        else
            usage();
    }
    try{
        GCLvectorfield3d image1(f1);
        double nrmi1=L2norm(image1);
        GCLvectorfield3d image2(f2);
        double nrmi2=L2norm(image2);
        image2*=(-1.0);
        image1+=image2;
        double nrmdiff=L2norm(image1);
        cout << "L2 norm of vectors in file "<<f1<<"="<<nrmi1<<endl;
        cout << "L2 norm of vectors in file "<<f2<<"="<<nrmi2<<endl;
        cout << "L2 norm of difference="<<nrmdiff<<endl;
        if(write_output_file)
            image1.save(ofile,string("."));
    }catch(std::exception& err)
    {
        cerr << err.what()<<endl;
    }
}
