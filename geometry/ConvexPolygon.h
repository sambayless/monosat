/*
 * ConvexPolygon.h
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#ifndef CONVEXPOLYGON_H_
#define CONVEXPOLYGON_H_
#include "Polygon.h"

/**
 * A concrete, convex polygon (or, for D>2, a polytope)
 */
template<unsigned int D,class T=double>
class ConvexPolygon:public Polygon<D,T>{
public:

	~ConvexPolygon(){};
	bool contains(Point<D,T> & point);
	bool intersects(Shape<D,T> & s){
		return false;
	}

};


template<>
bool ConvexPolygon<2,double>::contains(Point<2,double> & point);

#endif /* CONVEXPOLYGON_H_ */
