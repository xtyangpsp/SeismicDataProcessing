/* This program is a format converter to use the raw text file that defines
   ak135 (http://rses.anu.edu.au/seismology/ak135/ak135t.html) and convert
   it to the "mod1d" format used by pwmig.  
 
 This program is a unix filter.  Input is the raw text file for ak135
 (this includes a one line header) and the output is the mod1d format 
 used by pwmig*/
#include <iostream>
#include <string>
#include <math.h>
using namespace std;
int main(int argc,char **argv)
{
    double z,vp,vs,rho;
    double z0,vp0,vs0,rho0;
    double gradP,gradS,gradd,dz;
    string modname("ak135");
    string vpname("Pvelocity");
    string vsname("Svelocity");
    string dname("Density");
    const double mindz(0.01);  
    // first line of the file here has labels so skip it
    char line[80];
    /* Skip first line*/
    cin.getline(line,80);
    cin>>z0;  cin>>vp0;  cin>>vs0;   cin>>rho0;
    do
    {
        cin>>z;  cin>>vp;  cin>>vs;   cin>>rho;
        /* Here we always used the gradient fit between points.   
           When we hit a discontinuity we add an extra point 
           to define a 10 m thick discontinuity - a bit artibrary 
           but should work with pwmig.*/
        if(fabs(z-z0)<mindz)
        {
            gradP=(vp-vp0)/mindz;
            gradS=(vs-vs0)/mindz;
            gradd=(rho-rho0)/mindz;
            z+=mindz;
        }
        else
        {
            gradP=(vp-vp0)/(z-z0);
            gradS=(vs-vs0)/(z-z0);
            gradd=(rho-rho0)/(z-z0);
        }
        cout << modname<<" "<<vpname<<" "<<z0<<" "<<vp0<<" "<<gradP<<endl;
        cout << modname<<" "<<vsname<<" "<<z0<<" "<<vs0<<" "<<gradS<<endl;
        cout << modname<<" "<<dname<<" "<<z0<<" "<<rho0<<" "<<gradd<<endl;
        z0=z;
        vp0=vp;
        vs0=vs;
        rho0=rho;
    }while(!cin.eof());   // This will drop the last point, but irrelevant for pwmig
}
