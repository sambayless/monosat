/*
 * Polygon.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#include "ConvexPolygon.h"
#include "mtl/Vec.h"
template<>
bool ConvexPolygon<2,mpq_class>::containsInRange(const Point<2,mpq_class> & point, int firstVertex,int lastVertex){
	//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
	//this is correct _only_ for convex polygons
	 std::vector<Point<2,mpq_class> > &  w = getVertices();
	 dbg_orderClockwise();
	 if(w.size()==0)
		 return false;
	 else if(w.size()==1){
		 return w[0]==point;
	 }/*else if (w.size()==2){
		 //true if the point is between (inclusive) the other two points.

	 }*/
	 if(lastVertex<0)
		 lastVertex=w.size()-1;
	 //note: this can also compute the area (which is the sum of p2[0]*p1[1] - p1[0]*p2[1]); could potentially combine these...
	 for(int i = firstVertex;i!=lastVertex;i=(i+1% w.size() )){
		 Point<2,mpq_class> p1 = i>0 ? (w[i-1]-point):((w[lastVertex]-point));
		 Point<2,mpq_class> p2 = w[i]-point;
		 bool contained = (p2[0]*p1[1] - p1[0]*p2[1]) >=0;
		 if(!contained){
			 return false;
		 }
	 }
	 return true;
}



template<>
bool ConvexPolygon<2,mpq_class>::intersects(Shape<2,mpq_class> & shape){
	if(shape.getType()==CONVEX_POLYGON){
		ConvexPolygon & c = *(ConvexPolygon*) &shape;
		//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
		 std::vector<Point<2,mpq_class> > &  w = getVertices();
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
	assert(false);
	return false;
}
template<>
bool ConvexPolygon<2,mpq_class>::contains(const Point<2,mpq_class> & point){
	return containsInRange(point,0,this->size()-1);
}
template<>
bool ConvexPolygon<2,double>::containsInRange(const Point<2,double> & point, int firstVertex,int lastVertex){
	//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
	//this is correct _only_ for convex polygons
	 std::vector<Point<2,double> > &  w = getVertices();
	 dbg_orderClockwise();
	 if(w.size()==0)
		 return false;
	 else if(w.size()==1){
		 return w[0]==point;
	 }/*else if (w.size()==2){
		 //true if the point is between (inclusive) the other two points.

	 }*/
	 if(lastVertex<0)
		 lastVertex=w.size()-1;
	 //note: this can also compute the area (which is the sum of p2[0]*p1[1] - p1[0]*p2[1]); could potentially combine these...
	 for(int i = firstVertex;i!=lastVertex;i=(i+1% w.size() )){
		 Point<2,double> p1 = i>0 ? (w[i-1]-point):((w[lastVertex]-point));
		 Point<2,double> p2 = w[i]-point;
		 bool contained = (p2[0]*p1[1] - p1[0]*p2[1]) >=0;
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
		 std::vector<Point<2,double> > &  w = getVertices();
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
bool ConvexPolygon<1,double>::containsInRange(const Point<1,double> & point, int firstVertex,int lastVertex){
	return false;
}
template<>
bool ConvexPolygon<1,double>::contains(const Point<1,double> & point){
	return containsInRange(point,0,this->size()-1);
}
template<>
bool ConvexPolygon<2,double>::contains(const Point<2,double> & point){
	return containsInRange(point,0,this->size()-1);
}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::contains(const Point<D,T> & point)
{
	return containsInRange(point,0,this->size()-1);
}


template<>
bool ConvexPolygon<1,double>::intersects(Shape<1,double> & shape){
	return false;
}
