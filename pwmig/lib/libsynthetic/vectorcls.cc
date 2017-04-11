#include <iostream>
#include <cmath>
#include "vectorcls.h"

using namespace std;

vectorcls& vectorcls::operator+=(const vectorcls &lhs){
	v[0]+=lhs.v[0];
	v[1]+=lhs.v[1];
	v[2]+=lhs.v[2];
	return *this;
}
vectorcls& vectorcls::operator-=(const vectorcls &lhs){
	v[0]-=lhs.v[0];
	v[1]-=lhs.v[1];
	v[2]-=lhs.v[2];
	return *this;
}
//vectorcls& operator*=(const vectorcls &lhs);
	
vectorcls vectorcls::operator +(const vectorcls &lhs){
	
	return vectorcls(*this)+=lhs;
}
vectorcls vectorcls::operator -(const vectorcls &lhs){

	return vectorcls(*this)-=lhs;
}
vectorcls vectorcls::operator *(const vectorcls &lhs){
	//double v0, v1, v2;
	return vectorcls(v[1]*lhs.v[2]-v[2]*lhs.v[1],
		     v[2]*lhs.v[0]-v[0]*lhs.v[2],
		     v[0]*lhs.v[1]-v[1]*lhs.v[0]);
}

double dotproduct(const vectorcls &lhs, const vectorcls &rhs){
	return (rhs.v[0]*lhs.v[0]+rhs.v[1]*lhs.v[1]+rhs.v[2]*lhs.v[2]);
	
	
}
double modul(const vectorcls &lhs){
	return sqrt(dotproduct(lhs, lhs));
	return sqrt(lhs.v[0]*lhs.v[0]+lhs.v[1]*lhs.v[1]+lhs.v[2]*lhs.v[2]);
}
vectorcls& vectorcls::operator/=(const double &rhs){
	v[0]/=rhs;
	v[1]/=rhs;
	v[2]/=rhs;
	return *this;

}
vectorcls& vectorcls::operator*=(const double &rhs){
//times a scaler constant.
        v[0]*=rhs;
        v[1]*=rhs;
        v[2]*=rhs;
        return *this;

}
