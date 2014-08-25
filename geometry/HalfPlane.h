/*
 * Line.h
 *
 *  Created on: Jun 25, 2014
 *      Author: sam
 */
#ifndef HALF_PLANE_H_
#define HALF_PLANE_H_
#include <math.h>
#include "Shape.h"
#include "mtl/Vec.h"
#include "ConvexPolygon.h"
#include "Polygon.h"
#include "bounds/BoundingBox.h"
#include "Line.h"
#include "LineSegment.h"
/**
 * A half plane, defined by a (possibly not unit length) normal, and a distance along the normal from origin.
 */
template<unsigned int D,class T>
class HalfPlane:public Shape<D,T>{
public:
	Point<D,T> normal;
	T distance;

	bool isSingular(){
		return normal==0;
	}

	HalfPlane(){

	}
	HalfPlane(const Point<D,T> & normal, const T & distance):normal(normal),distance(distance){

	}

	ShapeType getType(){
		return HALF_PLANE;
	}

	bool contains(const Point<D,T> & point, bool inclusive){
		 T projection = normal.dot(point);
		 if(inclusive){
			 return projection<=distance;
		 }else{
			 return projection<distance;
		 }
	}
	bool intersects(Shape<D,T> & s, bool inclusive){
		if(s.getType()==LINE_SEGMENT){
			//a line segment intersects the half plane if its end points are on different sides of the half plane.
			LineSegment<D,T> & other = (LineSegment<D,T> & ) s;
			int side_a = this->whichSide(other.a);
			int side_b = this->whichSide(other.a);
			if(inclusive){
				return side_a==0 || side_b==0 || (side_a!=side_b);
			}else{
				return (side_a<0 && side_b>0) || (side_b<0 && side_a>0);
			}
		}else if(s.getType()==LINE){
			//Either the line is exactly parallel to the half plane, or it must intersect it eventually.
			Line<D,T> & other = (Line<D,T> & ) s;

			Point<D,T> other_dir = other.b - other.a;

			if(normal.dot( other_dir)==0){
				//if other dir is _exactly_ orthogonal to the normal of this plane, then the line is collinear to the half plane.
				//in that case, we need to check which side of the half plane the line is.
				//so pick _any_ point on the line and test it.
				int side = this->whichSide(other.a);
				if (side==0){
					return inclusive;
				}else{
					return side<0;
				}
			}else{
				//the line is not parallel to the half plane, and so must intersect eventually
				return true;
			}
		}else if(s.getType()==CONVEX_POLYGON || s.getType()==POLYGON){
			//a half plane intersects a polygon if any point of the convex polygon is inside the half plane.
			//can this be done more efficiently? In log time, possibly?
			Polygon<D,T> & poly = (Polygon<D,T> &)s;
			for(auto& p:poly){
				T side = whichSide(p);
				if(side==0 && inclusive){
					return true;
				}else if (side<0){
					return true;
				}
			}
			return false;
		}/*else if (s.getType()==BOUNDING_BOX){
			AbstractBoundingBox<D,T> & poly = (AbstractBoundingBox<2,T> &)s;
			for(auto& p:poly){
				T side = whichSide(p);
				if(side==0 && inclusive){
					return true;
				}else if (side<0){
					return true;
				}
			}
			return false;
		}*//*else if (s.getType()==BOUNDING_SPHERE){
			BoundingSphere<D,T> & sphere = (BoundingSphere<2,T> &)s;
			 T projection = normal.dot(sphere.center);
			 if(inclusive){
				 return projection-sphere.radius <=distance;
			 }else{
				 return projection-sphere.radius<distance;
			 }
		}*/
		assert(false);
	}

	//> 0 if the point is 'right' of the half plane, <0 if 'left' of the half plane, 0 if exactly on the half plane.
	int whichSide(const Point<D,T> & point){
		 T projection = normal.dot(point);
		 if(projection==distance){
			 return 0;
		 }else if (projection>distance){
			 return 1;
		 }else{
			 return -1;
		 }
	}
};

#endif



