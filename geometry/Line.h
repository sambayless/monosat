/*
 * Line.h
 *
 *  Created on: Jun 25, 2014
 *      Author: sam
 */
#ifndef LINE_H_
#define LINE_H_
#include <math.h>
#include "Shape.h"
#include "mtl/Vec.h"
#include "ConvexPolygon.h"
/**
 * A line
 */
template<unsigned int D,class T>
class Line:public Shape<D,T>{
public:
	//A line is defined by two non-equal points that it passes through
	Point<D,T> a;
	Point<D,T> b;

	bool isSingular(){
		return a==b;
	}

	Line(){

	}
	Line(const Point<D,T> & a, const Point<D,T> & b):a(a),b(b){

	}

	ShapeType getType(){
		return LINE;
	}

	bool contains(const Point<D,T> & point);
	bool intersects(Shape<D,T> & s);

	//> 0 if the point is 'right' of the line, <0 if 'left' of the line, 0 if exactly on the line.
	int whichSide(const Point<D,T> & point);

};

template<class T>
class Line<2,T>:public Shape<2,T>{
public:
	Point<2,T> a;
	Point<2,T> b;

	bool isSingular(){
		return a==b;
	}

	Line(){

	}

	Line(const Point<2,T> & a, const Point<2,T> & b):a(a),b(b){

	}

	ShapeType getType(){
		return LINE;
	}
	bool contains(const Point<2,T> & point);
	bool intersects(Shape<2,T> & s);

	bool intersects(Line<2,T> & other, Point<2,T> & intersection, bool & colinear);
	//> 0 if the point is 'right' of the line, <0 if 'left' of the line, 0 if exactly on the line.
	int whichSide(const Point<2,T> & point){
		T val= crossDif(b,a,point);
		if(val==0)
			return 0;
		return val>0?1:-1;
		//return ((b.x - a.x)*(point.y - a.y) - (b.y - a.y)*(point.x - a.x));
	}


};

template<class T>
bool Line<2,T>::contains(const Point<2,T> & point){
	return crossDif(a,b,point)==0;//is this correct?
}

template<class T>
bool Line<2,T>::intersects(Line<2,T> & other, Point<2,T> & intersection, bool & collinear){

	//from http://stackoverflow.com/a/565282
	collinear=false;
	bool intersecting=false;
	auto & p = a;
	Point<2,T> r = b-a;
	Point<2,T> & q = other.a;
	Point<2,T> s = other.b - other.a;

	T rs = cross2d(r,s);
	T qpr = cross2d(q-p,r);

	if(eq_epsilon(rs)){
		if(eq_epsilon(qpr)){
			//If r × s = 0 and (q − p) × r = 0, the two lines are collinear
			collinear=true;
		}else{
			//If r × s = 0 and (q − p) × r ≠ 0, then the two lines are parallel and non-intersecting.
		}
	}else{
		T t = cross2d((q - p),  s) / (cross2d(r , s));
		intersecting=true;
		intersection =r;
		intersection*=t;
		intersection+=p;
	}
	return collinear || intersecting;
}
template<class T>
bool Line<2,T>::intersects(Shape<2,T> & s){
	if(s.getType()==LINE){
		//Is this correct? This hasn't been tested.
		Line<2,T> & other = (Line<2,T> & ) s;
		T side1 = crossDif(a,b,other.a);
		if(side1==0)
			return true;//point is exactly on the line
		T side2 = crossDif(a,b,other.b);
		if(side2==0)
			return true;//point is exactly on the line
		if(side1!=side2){
			return true;//the lines must intersect;
		}
		return false;
	}else if(s.getType()==CONVEX_POLYGON){
		//a line intersects a convex polygon if it can find two points on opposite sides of the line.
		//there may be more efficient tests we could apply.
		Polygon<2,T> & poly = (Polygon<2,T> &)s;
		bool found_left = false;
		bool found_right = false;
		for(auto& p:poly){
			T side = whichSide(p);
			if(side==0){
				return true;
			}else if (side>0){
				found_right=true;
				if(found_left)
					return true;
			}else{
				found_left=true;
				if(found_right)
					return true;
			}
		}
		return false;
	}
	return false;
}

#endif



