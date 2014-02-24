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

using namespace Minisat;


template<unsigned int D, class T=double>
struct Point{
	int size()const{
		return D;
	}
	T vector[D];
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
    Point(){
    	for(int i = 0;i<D;i++){
    		new (&vector[i]) T();
    	}
    }

    //Copy constructor
    Point(const Point<D,T> & P){
    	for(int i = 0;i<D;i++){
    		new (&vector[i]) T(P[i]);
    	}
    }

    Point( std::initializer_list<T> list ){
    	assert(list.size()==size());
    	vector=list;
    }

    T dot(const Point<D,T> & other){
    	T sum=T();
    	for (int i = 0;i<D;i++){
    		sum += vector[i]*other[i];
    	}
    	return sum;
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
