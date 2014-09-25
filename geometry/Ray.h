/*
 * Line.h
 *
 *  Created on: Jun 25, 2014
 *      Author: sam
 */
#ifndef RAY_H_
#define RAY_H_
#include <cmath>
#include "Shape.h"
#include "mtl/Vec.h"
#include "Line.h"
#include <algorithm>
#include "ConvexPolygon.h"
/**
 * A ray defined by one end-point and a point it passes through
 */
template<unsigned int D,class T>
class Ray:public ConvexPolygon<D,T>{
public:
	//A ray segment is defined by an end-point and a point that it passes through
	Point<D,T> a;
	Point<D,T> b;



	Ray(){

	}
	Ray(const Point<D,T> & end_point, const Point<D,T> & b):a(end_point),b(b){

	}

	ShapeType getType(){
		return RAY;
	}

	bool contains(const Point<D,T> & point, bool inclusive);
	bool intersects(Shape<D,T> & s);

	//> 0 if the point is 'right' of the line, <0 if 'left' of the line, 0 if exactly on the line.
	int whichSide(const Point<D,T> & point);
	bool collinear(const  Point<D,T> & a,  const Point<D,T> &b){
		return contains(a, true) || contains(b,true);
	}
	int size()const{
		return 2;
	}
	void update(){

	}

	const Point<D,T>& operator [] (int index) const {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		if(index==0){
			return a;
		}else{
			assert(index==1);
			return b;
		}
	}
	Point<D,T>&       operator [] (int index)       {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		if(index==0){
			return a;
		}else{
			assert(index==1);
			return b;
		}
	}

};

template<class T>
class Ray<2,T>:public ConvexPolygon<2,T>{
public:
	Point<2,T> a;
	Point<2,T> b;


	Ray(){

	}

	Ray(const Point<2,T> & a, const Point<2,T> & b):a(a),b(b){

	}

	ShapeType getType(){
		return RAY;
	}
	bool contains(const Point<2,T> & point, bool inclusive);
	bool intersects(Shape<2,T> & s, bool inclusive);
	bool intersects(Shape<2,T> & shape, Point<2,T> & intersection, bool inclusive);
	bool intersects(LineSegment<2,T> & line, Point<2,T> & intersection, bool & overlapping, bool inclusive, bool calc_intersection);
	bool intersects(Line<2,T> & line, Point<2,T> & intersection, bool & overlapping, bool inclusive, bool calc_intersection);
	bool intersects(Ray<2,T> & other, Point<2,T> & intersection, bool & overlapping, bool inclusive, bool calc_intersection);
	//> 0 if the point is 'right' of the line, <0 if 'left' of the line, 0 if exactly on the line.
	int whichSide(const Point<2,T> & point){
		T val= crossDif(b,a,point);
		if(val==0)
			return 0;
		return val>0?1:-1;
		//return ((b.x - a.x)*(point.y - a.y) - (b.y - a.y)*(point.x - a.x));
	}
	bool collinear(const  Point<2,T> & a,  const Point<2,T> &b){
		return contains(a, true) || contains(b,true);
	}
	int size()const{
		return 2;
	}
	void update(){

	}

	const Point<2,T>& operator [] (int index) const {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		if(index==0){
			return a;
		}else{
			assert(index==1);
			return b;
		}
	}
	Point<2,T>&       operator [] (int index)       {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		if(index==0){
			return a;
		}else{
			assert(index==1);
			return b;
		}
	}

private:

};

template<class T>
bool Ray<2,T>::contains(const Point<2,T> & point, bool inclusive){
	if(!inclusive)
		return false;
	//adapted from http://stackoverflow.com/a/328122
	if(!eq_epsilon(crossDif(a,b,point))){
		return false;
	}
	//point is on the line; now check if it is in the ray.
	//it is in the ray if point.dot(b-a)>=0
	static Point<2,T> dif;
	dif = b;
	dif -=a;
	return point.dot(dif)>=0;
}

template<class T>
bool Ray<2,T>::intersects(Line<2,T> & other, Point<2,T> & intersection, bool & overlapping, bool inclusive, bool calc_intersection=true){

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
		if(0<=t){// && t<=1 ){//&& 0<= u && u<=1
			//If r × s ≠ 0 and 0 ≤ t ≤ 1 and 0 ≤ u ≤ 1, the two line segments meet at the point p + t r = q + u s.
			intersecting=true;
			if(calc_intersection){
				//intersection = p + (r*t);
				intersection =r;
				intersection*=t;
				intersection+=p;
			}
		}else{
			//Otherwise, the line segment and the ray are not parallel but do not intersect.
		}
	}
	return overlapping || intersecting;
}
template<class T>
bool Ray<2,T>::intersects(LineSegment<2,T> & other, Point<2,T> & intersection, bool & overlapping, bool inclusive, bool calc_intersection){

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
			//If r × s = 0 and (q − p) × r = 0, then the line segment and the ray are collinear.
			//If in addition, either 0 <= (q − p) · r <= r · r or 0 <= (p − q) · s ≤ s · s, then the two lines are overlapping.
			//How should overlapping line segments be treated for the purposes of inclusive/exclusive intersectiong?
			if(!inclusive)
				return false;

			//this test here needs to be modified for rays.
			//is this correct? its untested right now.

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
		T u = cross2d((q - p), r) / cross2d(r, s);

		// && t<=1 condition removed, because the ray extends to infinity.
		if(inclusive && 0<=t&& 0<= u && u<=1){// && t<=1
			//If r × s ≠ 0 and 0 ≤ t and 0 ≤ u ≤ 1, then the ray and line segment meet at the point p + t r = q + u s.
			intersecting=true;
			if(calc_intersection){
			//intersection = p + (r*t);
			intersection =r;
			intersection*=t;
			intersection+=p;
			}
		}else if(!inclusive && 0<t  && 0< u && u<1){// && t<1
			intersecting=true;
			if(calc_intersection){
			//intersection = p + (r*t);
			intersection =r;
			intersection*=t;
			intersection+=p;
			}
		}else{
			//Otherwise, the ray and line segment are not parallel but do not intersect.
		}
	}
	if(inclusive)
		return overlapping || intersecting;
	else
		return intersecting;
}

template<class T>
bool Ray<2,T>::intersects(Ray<2,T> & other, Point<2,T> & intersection, bool & overlapping, bool inclusive, bool calc_intersection){

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
			//If r × s = 0 and (q − p) × r = 0, then the line segment and the ray are collinear.
			//If in addition, either 0 <= (q − p) · r <= r · r or 0 <= (p − q) · s ≤ s · s, then the two lines are overlapping.
			//How should overlapping line segments be treated for the purposes of inclusive/exclusive intersectiong?
			if(!inclusive)
				return false;

			//this test here needs to be modified for rays.
			//is this correct? its untested right now.

			bool overlapping=false;
			T test = (q-p).dot(r);
			//this is utterly untested.
			if(0<=test ){//&& test<r.dot(r)
				overlapping=true;
			}else{
				//If r × s = 0 and (q − p) × r = 0, but neither 0 ≤ (q − p) · r ≤ r · r nor 0 ≤ (p − q) · s ≤ s · s, then the two lines are collinear but disjoint.
			}
		}else{
			//If r × s = 0 and (q − p) × r ≠ 0, then the two lines are parallel and non-intersecting.
		}
	}else{

		T t = cross2d((q - p),  s) / (cross2d(r , s));
		T u = cross2d((q - p), r) / cross2d(r, s);

		// && t<=1 condition removed, because the ray extends to infinity.
		if(inclusive && 0<=t&& 0<= u ){// && t<=1 && u<=1
			//If r × s ≠ 0 and 0 ≤ t and 0 ≤ u ≤ 1, then the ray and line segment meet at the point p + t r = q + u s.
			intersecting=true;
			if(calc_intersection){
			//intersection = p + (r*t);
			intersection =r;
			intersection*=t;
			intersection+=p;
			}
		}else if(!inclusive && 0<t  && 0< u){// && t<1  && u<1
			intersecting=true;
			if(calc_intersection){
			//intersection = p + (r*t);
			intersection =r;
			intersection*=t;
			intersection+=p;
			}
		}else{
			//Otherwise, the rays are not parallel but do not intersect.
		}
	}
	if(inclusive)
		return overlapping || intersecting;
	else
		return intersecting;
}
template<class T>
bool Ray<2,T>::intersects(Shape<2,T> & s, bool inclusive){
	if(s.getType()==RAY){
			bool ignore_b=false;
			static Point<2,T> ignore_point;
			return intersects(Ray<2,T> (s),ignore_point,ignore_b,inclusive,false);
	}else if(s.getType()==LINE){
		bool ignore_b=false;
		static Point<2,T> ignore_point;
		return intersects(Line<2,T> (s),ignore_point,ignore_b,inclusive,false);
	}else	if(s.getType()==LINE_SEGMENT){
		bool ignore_b=false;
		static Point<2,T> ignore_point;
		return intersects(LineSegment<2,T> (s),ignore_point,ignore_b,inclusive,false);
	}else if(s.getType()==CONVEX_POLYGON){
		//see http://stackoverflow.com/a/312332 for more efficient implementations.

		//a line intersects a convex polygon if it can find two points on opposite sides of the line.
		//there may be more efficient tests we could apply.
		Polygon<2,T> & poly = (Polygon<2,T> &)s;
		bool found_left = false;
		bool found_right = false;
		for(auto& p:poly){
			T side = whichSide(p);
			if(side==0 && inclusive){
				return true;
			}else if (side>0){
				found_right=true;
				if(found_left)
					return true;
			}else if (side<0){
				found_left=true;
				if(found_right)
					return true;
			}
		}
		return false;
	}
	return false;
}

template<class T>
bool Ray<2,T>::intersects(Shape<2,T> & shape, Point<2,T> & intersection, bool inclusive){
	if(shape.getType()==LINE){
		bool ignore;
		return this->intersects((Line<2,T> &)shape,intersection,ignore,inclusive);
	}else if(shape.getType()==LINE_SEGMENT){
		bool ignore;
		return this->intersects((Line<2,T> &)shape,intersection,ignore,inclusive);
	}else if(shape.getType()==RAY){
		bool ignore;
		return this->intersects((Ray<2,T> &)shape,intersection,ignore,inclusive);
	}else if (shape.getType()==CONVEX_POLYGON ||shape.getType()==POLYGON){
		Polygon<2,T> & polygon = (Polygon<2,T> &)shape;

		//find the point closest to the origin of the ray that is in the polygon (if any).
		//first check if the origin itself is contained
		if(polygon.contains(a,inclusive)){
			intersection=a;
			return true;
		}

		//ok, now test each line segment of the polygon against the ray
		static LineSegment<2,T> check;
		static Point<2,T> point;
		static Point<2,T> line_test;
		bool found=false;
		line_test = b-a;
		static T least_distance= numeric<T>::inf;
		for (int i = 0;i<polygon.size();i++){
			check.a = polygon[i-1];
			check.b = polygon[i];
			if(intersects(check,point,inclusive)){
				T distance = line_test.dot(point);
				assert(distance>=0);
				if(distance<least_distance){
					least_distance=distance;
					intersection=point;
					found=true;
				}
			}
		}
		return found;
	}
}


template<unsigned int D,class T>
std::ostream & operator<<(std::ostream & str, Ray<D,T> const & p){
	str << "Ray" << p.a <<","<<p.b<<"]";
	return str;
}
#endif



