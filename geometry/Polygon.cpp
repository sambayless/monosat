/*
 * Polygon.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#include "Polygon.h"
#include "mtl/Vec.h"

//Note that this is subject to rounding errors. It also might only be correct for convex polygons.
//Also, this is only correct for planar (non-self-intersecting) polygons.
template<>
double Polygon<2,double>::getArea(){
	 vec<Point2D> &  w = getVertices();
	double sum = 0;
	for (int i = 1;i<w.size();i++){
		Point2D prev = w[i-1];
		Point2D cur = w[i];
		sum += prev[0]*cur[1]-cur[0]*prev[1];
	}
	return sum/2.0;
}
//Note that this is subject to rounding errors.
template<>
double Polygon<2,double>::getPerimeter(){
	vec<Point2D> & w = getVertices();
	double sum = 0;
	for (int i = 1;i<w.size();i++){
		Point2D prev = w[i-1];
		Point2D cur = w[i];
		double xdist = cur[0]-prev[0];
		double ydist=cur[1]-prev[1];
		sum += sqrt(xdist*xdist + ydist*ydist);
	}
	return sum;
}
