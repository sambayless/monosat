/*
 * GraphTheoryTypes.h
 *
 *  Created on: 2014-01-08
 *      Author: sam
 */

#ifndef GEOMETRY_TYPES_H_
#define GEOMETRY_TYPES_H_
#include <initializer_list>
#include "core/SolverTypes.h"
#include <cmath>
#include <algorithm>
using namespace Minisat;


template<unsigned int D, class T=double>
struct Point{
	int size()const{
		return D;
	}
	T vector[D];
	T & x;
	T & y;
	T & z;

    // Vector interface:
    const T& operator [] (int index) const {assert(index<D);assert(index>=0); return vector[index];}
    T&       operator [] (int index)       {assert(index<D);assert(index>=0);  return vector[index];}

/*    bool operator == (const Point<D,T> & other)const{
    	for(int i = 0;i<D;i++){
    		if(vector[i]!= other[i])
    			return false;
    	}
    	return true;
    }*/
    Point():x(vector[0]),y(vector[1]),z(vector[2]){
    	for(int i = 0;i<D;i++){
    		new (&vector[i]) T();
    	}
    }
    Point(const vec<T> & list  ):x(vector[0]),y(vector[1]),z(vector[2]){
    	for(int i = 0;i<D;i++){
    		vector[i] = list[i];
    	}
    }
    //Copy constructor
    Point(const Point<D,T> & v):x(vector[0]),y(vector[1]),z(vector[2]){
    	assert(v.size()==D);
    	for(int i = 0;i<D;i++){
    		new (&vector[i]) T(v[i]);
    	}
    }

    Point( std::initializer_list<T> list ):x(vector[0]),y(vector[1]),z(vector[2]){
    	assert(list.size()==size());
    	vector=list;
    }
    template<typename... Ts>
    Point( Ts... args ):vector{args...},x(vector[0]),y(vector[1]),z(vector[2]){

    }
    Point& operator=(const Point & v)
    {
    	for(int i = 0;i<D;i++){
			new (&vector[i]) T(v[i]);
		}
      return *this;
    }

    T dot(const Point<D,T> & other){
    	T sum=T();
    	for (int i = 0;i<D;i++){
    		sum += vector[i]*other[i];
    	}
    	return sum;
    }
    void zero(){
    	for(int i = 0;i<D;i++){
    		vector[i]=0;
    	}
    }
    Point<D,T> &       operator += (const Point<D,T>& other) {
    	for(int i = 0;i<D;i++){
    		vector[i]+=other[i];
    	}
    	return *this;
    }
    Point<D,T> &       operator -= (const Point<D,T>& other) {
    	for(int i = 0;i<D;i++){
    		vector[i]-=other[i];
    	}
    	return *this;
    }
    Point<D,T> &       operator /= (const T & scalar) {
    	for(int i = 0;i<D;i++){
    		vector[i]/=scalar;
    	}
    	return *this;
    }
    Point<D,T> &       operator *= (const T & scalar) {
    	for(int i = 0;i<D;i++){
    		vector[i]*=scalar;
    	}
    	return *this;
    }
    Point<D,T>        operator * (const T & scalar)const {
    	Point<D,T> ret;
    	for(int i = 0;i<D;i++){
    		ret[i]= vector[i]* scalar;
    	}
    	return ret;
    }
    Point<D,T>        operator / (const T & scalar)const {
    	Point<D,T> ret;
    	for(int i = 0;i<D;i++){
    		ret[i]= vector[i]/ scalar;
    	}
    	return ret;
    }
    //Note that this computes a DOUBLE. Should enforce that is a safe underapproximation of the real distance between these points (that is, the actual distance is guaranteed to be >= the returned value)
    double distance_underapprox(Point<D,T> other){
    	T sum=T(0);
    	for(int i =0;i<D;i++){
    		T difference = (vector[i]-other[i]);
    		sum+= difference*difference;
    	}
    	double distance = sqrt((double)sum);
    	return distance;
    }
};
template<unsigned int D, class T=double>
inline bool operator==(const Point<D,T>& lhs, const Point<D,T>& rhs){
	for(int i = 0;i<D;i++){
		if(lhs[i]!= rhs[i])
			return false;
	}
	return true;
}

template<unsigned int D, class T=double>
inline Point<D,T> operator+(const Point<D,T> &a, const Point<D,T> &b)
{
	Point<D,T> p;
	for(int i = 0;i<D;i++){
		p[i]=a[i]+b[i];
	}
    return p;
}

template<unsigned int D, class T=double>
inline Point<D,T> operator-(const Point<D,T> &a, const Point<D,T> &b)
{
	Point<D,T> p;
	for(int i = 0;i<D;i++){
		p[i]=a[i]-b[i];
	}
    return p;
}



typedef Point<2,double> Point2D;
typedef Point<3,double> Point3D;

template<unsigned int D,class T=double>
struct SortBy{
	int sortOn;
	bool operator()(Point<D,T> & a,Point<D,T> & b)const{
		return a[sortOn]<b[sortOn];
	}
	SortBy(int dimensionToSort):sortOn(dimensionToSort){}
};

template<unsigned int D,class T=double>
struct SortLexicographic{

	bool operator()(Point<D,T> & a,Point<D,T> & b)const{
		for(int i =0;i<D;i++)
			if(a[i]<b[i])
				return true;
			else if (a[i]>b[i])
				return false;
		return false;
	}

};
    // Returns a random float 0 <= x < 1. Seed must never be 0.
    static inline double drand(double& seed) {
    assert(seed!=0);
       seed *= 1389796;
       int q = (int)(seed / 2147483647);
       seed -= (double)q * 2147483647;
       return seed / 2147483647; }

    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline int irand(double& seed, int size) {
       return (int)(drand(seed) * size); }




#endif /* GEOMETRY_TYPES_H_ */
