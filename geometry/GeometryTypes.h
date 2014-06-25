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
#include "mtl/Rnd.h"
#include <cmath>
#include <algorithm>

using namespace Minisat;

/**
 * c++14 version:
 * template <typename T> constexpr T epsilon;

template <> constexpr float  epsilon<float>  = 0.001f;
template <> constexpr double epsilon<double> = 0.000001;
 */
template <typename T> struct epsilon;

template <> struct epsilon<float>
{
private:
    float const value = 0.000001f;
public:

    operator double()const{return value;}
};

template <> struct epsilon<double>
{
private:
    double const value = 0.000000001;
public:

    operator double()const{return value;}
};




template<class T> bool equal_epsilon(T a, T b);

template<> inline bool equal_epsilon(double a, double b){
	return std::abs(a-b)<=epsilon<double>();
}
template<> inline bool equal_epsilon(float a, float b){
	return std::abs(a-b)<=epsilon<float>();
}

template<unsigned int D, class T=double>
struct Point{
	int size()const{
		return D;
	}
	int id=-1;//optional identifier variable
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
    Point():id(-1),x(vector[0]),y(vector[1]),z(vector[2]){
    	for(int i = 0;i<D;i++){
    		new (&vector[i]) T();
    	}
    }
    Point(const vec<T> & list  ):id(-1),x(vector[0]),y(vector[1]),z(vector[2]){
    	for(int i = 0;i<D;i++){
    		vector[i] = list[i];
    	}
    }
    //Copy constructor
    Point(const Point<D,T> & v):id(v.id),x(vector[0]),y(vector[1]),z(vector[2]){
    	assert(v.size()==D);
    	for(int i = 0;i<D;i++){
    		new (&vector[i]) T(v[i]);
    	}
    }

    Point( std::initializer_list<T> list ):id(-1),x(vector[0]),y(vector[1]),z(vector[2]){
    	assert(list.size()==size());
    	vector=list;
    }
    template<typename... Ts>
    Point( Ts... args ):id(-1),vector{args...},x(vector[0]),y(vector[1]),z(vector[2]){
    	int a =1;
    }
    Point& operator=(const Point & v)
    {
    	id=v.id;
    	for(int i = 0;i<D;i++){
			new (&vector[i]) T(v[i]);
		}
      return *this;
    }
    int getID()const{
    	assert(hasID());
    	return id;
    }
    bool hasID()const{
    	return id>=0;
    }
    void setID(int id){
    	this->id=id;
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


template<unsigned int D, class T>
class GeometryTheorySolver;


#endif /* GEOMETRY_TYPES_H_ */
