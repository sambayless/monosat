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
	Point<D,T> circleCenter;//should generalize this to an arbitrary bounding volume...
	T circleRadius;
	virtual ~Polygon(){};
	virtual ShapeType getType(){
		return POLYGON;
	}
	virtual bool contains(Point<D,T> & point){
		return false;
	}
	virtual bool intersects(Shape<D,T> & s){
		return false;
	}


	void update(){
		updateCircleBound();
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
	void updateCircleBound();
	virtual T getArea();
	virtual T getPerimeter();
};
template<>
double Polygon<2,double>::getArea();

template<>
double Polygon<2,double>::getPerimeter();

template<unsigned int D,class T>
void Polygon<D,T>::updateCircleBound(){
	circleCenter.zero();
	for(int i = 0;i<vertices.size();i++){
		circleCenter+=vertices[i];
	}
	circleCenter/=T(vertices.size());
	circleRadius=T(0);
	for(int i = 0;i<vertices.size();i++){
		T dist = circleCenter.distance( vertices[i]);
		if(dist>circleRadius){
			circleRadius=dist;
		}
	}
}

#endif /* POLYGON_H_ */
