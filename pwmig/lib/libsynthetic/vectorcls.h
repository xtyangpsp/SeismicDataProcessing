#ifndef _VECTORCLS_H_ 
#define	_VECTORCLS_H_ 
#include <iostream>

using namespace std;

class vectorcls{

public:
	vectorcls(){
	}
	vectorcls(double v0, double v1, double v2){
		v[0]=v0;
		v[1]=v1;
		v[2]=v2;
	}
	vectorcls(const vectorcls &cp){
		v[0]=cp.v[0];
		v[1]=cp.v[1];
		v[2]=cp.v[2];
	
	}
	
	vectorcls& operator=(const vectorcls &lhs){
		v[0]=lhs.v[0];
		v[1]=lhs.v[1];
		v[2]=lhs.v[2];
		return *this;
	}
	
	vectorcls& operator+=(const vectorcls &lhs);
	vectorcls& operator-=(const vectorcls &lhs);
	vectorcls& operator/=(const double &rhs);
	vectorcls& operator*=(const double &rhs);
	//vectorcls& operator*=(const vectorcls &lhs);
	
	vectorcls operator +(const vectorcls &lhs);
	vectorcls operator -(const vectorcls &lhs);
	vectorcls operator *(const vectorcls &lhs);
	
	friend double dotproduct(const vectorcls &lhs, const vectorcls &rhs);
	friend double modul(const vectorcls &lhs);
	//modul(const vectorcls) compute the modulus of a vector
	
	double* get_address (int startind){
		return &v[startind];
	}
        double* get_address (){
	// this is the default entry with no parameter required
                return &v[0];
        }

	friend ostream& operator<<(ostream& out, const vectorcls& vec) // output
	{
	    out << "(" << vec.v[0] << ", " << vec.v[1] <<", "<<vec.v[2] << ")";
	    return out;
	}
	friend istream& operator>>(istream& in, vectorcls& vec) // input
	{
		//double x, y;
		in >> vec.v[0] >> vec.v[1]>>vec.v[2];
		return in;
	}


private:

	double v[3];
	//a 3-component vector





};

#endif
