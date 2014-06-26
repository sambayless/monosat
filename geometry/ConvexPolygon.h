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

	virtual ~ConvexPolygon(){};
	virtual ShapeType getType(){
		return CONVEX_POLYGON;
	}
	bool contains(const Point<D,T> & point);
	bool containsInRange(const Point<D,T> & point, int firstVertex=0,int lastVertex=-1);
	bool intersects(Shape<D,T> & s);
};


template<>
bool ConvexPolygon<2,double>::containsInRange(const Point<2,double> & point, int firstVertex,int lastVertex);

template<>
bool ConvexPolygon<2,double>::intersects(Shape<2,double> & s);

template<>
bool ConvexPolygon<1,double>::containsInRange(const Point<1,double> & point, int firstVertex,int lastVertex);

template<>
bool ConvexPolygon<1,double>::intersects(Shape<1,double> & s);



#endif /* CONVEXPOLYGON_H_ */
