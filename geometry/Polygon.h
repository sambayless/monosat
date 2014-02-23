/*
 * Polygon.h
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#ifndef POLYGON_H_
#define POLYGON_H_
#include <math.h>
#include "Shape.h"

/**
 * A concrete polygon (or, for D>2, a polytope)
 */
template<unsigned int D,class T=double>
class Polygon:public Shape<D>{
public:
	//List of vertices in clockwise order
	vec<Point<D,T>> vertices;
	~Polygon(){};
	bool contains(Point<D,T> & point){
		return false;
	}
	bool intersects(Shape<D,T> & s){
		return false;
	}

	void clear(){
		vertices.clear();
	}
	int size(){
		return vertices.size();
	}
	//Returns the vertices of the polygon, in clockwise order.
	//Note: for convenience, the first point of the wrap is also the last point (it is duplicated).
	vec<Point<D,T> > & getVertices(){
		return vertices;
	}

	T getArea();
	T getPerimeter();
};
template<>
double Polygon<2,double>::getArea();

template<>
double Polygon<2,double>::getPerimeter();

#endif /* POLYGON_H_ */
