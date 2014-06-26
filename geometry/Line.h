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
/**
 * A line
 */
template<unsigned int D,class T=double>
class Line:public Shape<D>{
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

	bool contains(Point<D,T> & point);
	bool intersects(Shape<D,T> & s);

	//> 0 if the point is 'right' of the line, <0 if 'left' of the line, 0 if exactly on the line.
	T whichSide(Point<D,T> & point);

};

template<class T>
class Line<2,T>{
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
	bool contains(Point<2,T> & point);
	bool intersects(Shape<2,T> & s);

	//> 0 if the point is 'right' of the line, <0 if 'left' of the line, 0 if exactly on the line.
	T whichSide(Point<2,T> & point){
		return cross(b,a,point);
		//return ((b.x - a.x)*(point.y - a.y) - (b.y - a.y)*(point.x - a.x));
	}
private:
	T cross(const Point<2,T> &O, const Point<2,T>  &A, const Point<2,T>  &B)
	{
		return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
	}
};

template<class T>
bool Line<2,T>::contains(Point<2,T> & point){
	return cross(a,b,point)==0;//is this correct?
}

template<class T>
bool Line<2,T>::intersects(Shape<2,T> & s){
	if(s.getType()==LINE){
		//they intersect if not collinear

	}else if(s.getType()==CONVEX_POLYGON){
		//a line intersects a convex polygon if it can find two points on opposite sides of the line.
		//there may be more efficient tests we could apply.
		ConvexPolygon<2,T> & poly = (ConvexPolygon<2,T> &)s;
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
}

#endif



