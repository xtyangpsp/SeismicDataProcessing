#include <string>
#include <list>
#include <fstream>
#include <iostream>
using namespace std;
string extract_sta_name(char *filename)
{
    string path(filename);
    string sta;
    size_t right;
    right=path.rfind("/");
    if(right==string::npos)
        sta=path;
    else
        sta.assign(path,right+1,string::npos);
    return(sta);
}
void usage()
{
    cout << "MergeEarsModels < filenamelist > out"<<endl
        << "   filenamelist is list of model file names to merge to out"<<endl;
    exit(-1);
}
int main(int argc, char **argv)
{
    if(argc!=1) usage();
    char fname[256];
    while(cin.getline(fname,256))
    {
        string sta=extract_sta_name(fname);
        ifstream din;
        din.open(fname,ios::in);
        if(din.fail())
        {
            cerr << "open failed on filename="<<fname<<endl
               << "Will not continue.  Fix and rerun"<<endl;
            exit(-1);
        }
        char inputline[256];
        list<string> modeldata;
        modeldata.clear();
        while(din.getline(inputline,256))
        {
            modeldata.push_back(string(inputline));
        }
        cout << sta <<" "<<modeldata.size()<<endl;
        list<string>::iterator mptr;
        for(mptr=modeldata.begin();mptr!=modeldata.end();++mptr)
            cout << *mptr<<endl;
        din.close();
    }
}
