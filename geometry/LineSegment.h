/*
 * Line.h
 *
 *  Created on: Jun 25, 2014
 *      Author: sam
 */
#ifndef LINE_SEGMENT_H_
#define LINE_SEGMENT_H_
#include <cmath>
#include "Shape.h"
#include "mtl/Vec.h"
#include "Line.h"
#include <algorithm>
/**
 * A line segment defined by two end-points
 */
template<unsigned int D,class T>
class LineSegment:public Shape<D,T>{
public:
	//A line segment is defined by two end-points that it passes through
	Point<D,T> a;
	Point<D,T> b;



	LineSegment(){

	}
	LineSegment(const Point<D,T> & a, const Point<D,T> & b):a(a),b(b){

	}

	ShapeType getType(){
		return LINE_SEGMENT;
	}

	bool contains(const Point<D,T> & point, bool inclusive);
	bool intersects(Shape<D,T> & s);

	//> 0 if the point is 'right' of the line, <0 if 'left' of the line, 0 if exactly on the line.
	int whichSide(const Point<D,T> & point);

};

template<class T>
class LineSegment<2,T>:public Shape<2,T>{
public:
	Point<2,T> a;
	Point<2,T> b;


	LineSegment(){

	}

	LineSegment(const Point<2,T> & a, const Point<2,T> & b):a(a),b(b){

	}

	ShapeType getType(){
		return LINE_SEGMENT;
	}
	bool contains(const Point<2,T> & point, bool inclusive);
	bool intersects(Shape<2,T> & s, bool inclusive);

	bool intersects(LineSegment<2,T> & line, Point<2,T> & intersection, bool & overlapping, bool inclusive);
	bool intersects(Line<2,T> & line, Point<2,T> & intersection, bool & overlapping, bool inclusive);
	//> 0 if the point is 'right' of the line, <0 if 'left' of the line, 0 if exactly on the line.
	int whichSide(const Point<2,T> & point){
		T val= cross(b,a,point);
		if(val==0)
			return 0;
		return val>0?1:-1;
		//return ((b.x - a.x)*(point.y - a.y) - (b.y - a.y)*(point.x - a.x));
	}
private:

	static bool mightIntersect(LineSegment<2,T> & l1,LineSegment<2,T> &l2, bool inclusive);
};

template<class T>
bool LineSegment<2,T>::contains(const Point<2,T> & point, bool inclusive){
	if(!inclusive)
		return false;
	return false;
}

template<class T>
bool LineSegment<2,T>::intersects(Line<2,T> & other, Point<2,T> & intersection, bool & overlapping, bool inclusive){

	//from http://stackoverflow.com/a/565282
	overlapping=false;
	bool intersecting=false;
	auto & p = a;
	Point<2,T> r = b-a;
	Point<2,T> & q = other.a;
	Point<2,T> s = other.b - other.a;

	T rs = cross2d(r,s);
	T qpr = cross2d(q-p,r);


	if(eq_epsilon(rs)){
		if(eq_epsilon(qpr)){
			//If r × s = 0 and (q − p) × r = 0, then the two lines are collinear.
			//If in addition, either 0 <= (q − p) · r <= r · r or 0 <= (p − q) · s ≤ s · s, then the two lines are overlapping.
			bool overlapping=false;
			T test = (q-p).dot(r);
			if(0<=test && test<r.dot(r)){
				overlapping=true;
			}else{
				//If r × s = 0 and (q − p) × r = 0, but neither 0 ≤ (q − p) · r ≤ r · r nor 0 ≤ (p − q) · s ≤ s · s, then the two lines are collinear but disjoint.
			}
		}else{
			//If r × s = 0 and (q − p) × r ≠ 0, then the two lines are parallel and non-intersecting.
		}
	}else{

		T t = cross2d((q - p),  s) / (cross2d(r , s));
		//T u = cross2d((q - p), r) / cross2d(r, s);
		if(0<=t && t<=1 ){//&& 0<= u && u<=1
			//If r × s ≠ 0 and 0 ≤ t ≤ 1 and 0 ≤ u ≤ 1, the two line segments meet at the point p + t r = q + u s.
			intersecting=true;
			//intersection = p + (r*t);
			intersection =r;
			intersection*=t;
			intersection+=p;
		}else{
			//Otherwise, the two line segments are not parallel but do not intersect.
		}
	}
	return overlapping || intersecting;
}
template<class T>
bool LineSegment<2,T>::intersects(LineSegment<2,T> & other, Point<2,T> & intersection, bool & overlapping, bool inclusive){

	//from http://stackoverflow.com/a/565282
	overlapping=false;
	bool intersecting=false;
	auto & p = a;
	Point<2,T> r = b-a;
	Point<2,T> & q = other.a;
	Point<2,T> s = other.b - other.a;

	T rs = cross2d(r,s);
	T qpr = cross2d(q-p,r);


	if(eq_epsilon(rs)){
		if(eq_epsilon(qpr)){
			//If r × s = 0 and (q − p) × r = 0, then the two lines are collinear.
			//If in addition, either 0 <= (q − p) · r <= r · r or 0 <= (p − q) · s ≤ s · s, then the two lines are overlapping.
			bool overlapping=false;
			T test = (q-p).dot(r);
			if(0<=test && test<r.dot(r)){
				overlapping=true;
			}else{
				test = (p-q).dot(s);
				if(0<=test && test<= s.dot(s)){
					overlapping=true;
				}else{
					//If r × s = 0 and (q − p) × r = 0, but neither 0 ≤ (q − p) · r ≤ r · r nor 0 ≤ (p − q) · s ≤ s · s, then the two lines are collinear but disjoint.
				}
			}
		}else{
			//If r × s = 0 and (q − p) × r ≠ 0, then the two lines are parallel and non-intersecting.
		}
	}else{

		T t = cross2d((q - p),  s) / (cross2d(r , s));
		T u = cross2d((q - p), r) / cross2d(r, s);
		if(0<=t && t<=1 && 0<= u && u<=1){
			//If r × s ≠ 0 and 0 ≤ t ≤ 1 and 0 ≤ u ≤ 1, the two line segments meet at the point p + t r = q + u s.
			intersecting=true;
			//intersection = p + (r*t);
			intersection =r;
			intersection*=t;
			intersection+=p;
		}else{
			//Otherwise, the two line segments are not parallel but do not intersect.
		}
	}
	return overlapping || intersecting;
}

template<class T>
bool LineSegment<2,T>::mightIntersect(LineSegment<2,T> & l1,LineSegment<2,T> &l2, bool inclusive){

	T side1 = crossDif(l1.a,l1.b,l2.a);
	if(side1==0)
		return true;//point is exactly on the line
	T side2 = crossDif(l1.a,l1.b,l2.b);
	if(side2==0)
		return true;//point is exactly on the line
	if((side1 > 0) !=  (side2 > 0)){
		return true;//the lines might intersect;
	}
	Point<2,T> diff = l1.b-l1.a;
	Point<2,T> check = l2.a - l1.a;
	bool isARight = cross2d(diff,check)<0;

	check = l2.b - l1.a;
	bool isBRight = cross2d(diff,check)<0;

	if( isARight != isBRight){
		return true;
	}
	return false;
}

template<class T>
bool LineSegment<2,T>::intersects(Shape<2,T> & shape, bool inclusive){
	if(shape.getType()==LINE_SEGMENT){

		LineSegment<2,T> & other = (LineSegment<2,T> &)shape;
		//from http://martin-thoma.com/how-to-check-if-two-line-segments-intersect/

		//first, check whether bounding boxes intersect. This appears to be necessary for the correctness of this method (as opposed to just being an optimization)
		T minX = std::min(a.x, b.x);
		T minY = std::min(a.y,b.y);
		T maxX = std::max(a.x,b.x);
		T maxY = std::max(a.y,b.y);
		if(other.a.x < minX && other.b.x < minX){
			return false;
		}
		if(other.a.x > maxX && other.b.x > maxX){
			return false;
		}
		if(other.a.y < minY && other.b.y < minY){
			return false;
		}
		if(other.a.y > maxY && other.b.y > maxY){
			return false;
		}

		return mightIntersect(*this,other,inclusive) && mightIntersect(other,*this,inclusive);
	}else if(shape.getType()==CONVEX_POLYGON){
		return shape.intersects(*this, inclusive);
	}
	return false;
}
template<unsigned int D,class T>
std::ostream & operator<<(std::ostream & str, LineSegment<D,T> const & p){
	str << "Line=[" << p.a <<","<<p.b<<"]";
	return str;
}
#endif



