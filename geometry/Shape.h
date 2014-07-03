/*
 * Shape.h
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#ifndef SHAPE_H_
#define SHAPE_H_
#include <vector>
#include "GeometryTypes.h"

enum ShapeType{
	CONVEX_POLYGON,POLYGON,SHAPE,PLANE,LINE,LINE_SEGMENT
};
/**
 * A concrete shape
 */
template<unsigned int D,class T>
class Shape{

public:

	Shape(){
	}
	virtual ~Shape(){};
	int dimension(){
		return D;
	}

	virtual ShapeType getType(){
		return SHAPE;
	}

	virtual bool contains(const Point<D,T> & point)=0;
	virtual bool intersects(Shape<D,T> & s)=0;


};



#endif /* SHAPE_H_ */
