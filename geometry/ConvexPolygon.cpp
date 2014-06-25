/*
 * Polygon.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#include "ConvexPolygon.h"
#include "mtl/Vec.h"

template<>
bool ConvexPolygon<2,double>::contains(Point<2,double> & point){
	//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
	 vec<Point<2,double> > &  w = getVertices();
	 if(w.size()<3)
		 return false;
	 //note: this can also compute the area (which is the sum of p2[0]*p1[1] - p1[0]*p2[1]); could potentially combine these...
	 for(int i = 1;i<w.size();i++){
		 Point<2,double> p1 = w[i-1]-point;
		 Point<2,double> p2 = w[i]-point;
		 bool contained = (p2[0]*p1[1] - p1[0]*p2[1]) >0;
		 if(!contained){
			 return false;
		 }
	 }
	 return true;
}


template<>
bool ConvexPolygon<2,double>::intersects(Shape<2,double> & shape){
	if(shape.getType()==CONVEX_POLYGON){
		ConvexPolygon & c = *(ConvexPolygon*) &shape;
		//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
		 vec<Point<2,double> > &  w = getVertices();
		 //note: this can also compute the area (which is the sum of p2[0]*p1[1] - p1[0]*p2[1]); could potentially combine these...
	/*	 for(int i = 1;i<w.size();i++){
			 Point<2,double> p1 = w[i-1]-point;
			 Point<2,double> p2 = w[i]-point;
			 bool contained = (p2[0]*p1[1] - p1[0]*p2[1]) >0;
			 if(!contained){
				 return false;
			 }
		 }*/
	}
	 return true;
}
template<>
bool ConvexPolygon<1,double>::contains(Point<1,double> & point){
	return false;
}


template<>
bool ConvexPolygon<1,double>::intersects(Shape<1,double> & shape){
	return false;
}
