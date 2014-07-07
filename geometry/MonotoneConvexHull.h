/*
 * SimpleConvexHull.h
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */


#ifndef MONOTONE_CONVEXHULL_H_
#define MONOTONE_CONVEXHULL_H_
#include "ConvexHull.h"
#include "PointSet.h"
#include <gmpxx.h>
#include "bounds/BoundingBox.h"
#include <algorithm>
#include <iostream>

template<unsigned int D,class T>
class MonotoneConvexHull:public ConvexHull<D,T>{

	PointSet<D,T> & pointSet;
	ConvexPolygon<D,T> hull;
	long lastModification=-1;
public:
	MonotoneConvexHull(PointSet<D,T> & p):pointSet(p){
		hull.setBoundingVolume(new BoundingBox<D,T,Polygon<D,T>>(hull));
	}

	void update(){

		if(pointSet.getModifications()<=lastModification){
			return;
		}
		lastModification=pointSet.getModifications();
		if(D==2){
			update2d();
			return;
		}
		assert(false);
	}


	ConvexPolygon<D,T> & getHull(){
		update();
		return hull;
	}

private:
	T cross(const Point<D,T> &O, const Point<D,T> &A, const Point<D,T> &B);

	void update2d();
};

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
//from http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
template<>
inline double MonotoneConvexHull<2,double>::cross(const Point2D &O, const Point2D &A, const Point2D &B)
{
	return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
}


template<unsigned int D,class T>
void MonotoneConvexHull<D,T>::update2d(){


	std::vector<Point<D,T>> _points;
	pointSet.getEnabledPoints(_points);
	std::vector<Point<2,T>> & points = (std::vector<Point<2,T>> &)_points;

	std::sort(points.begin(),points.end(),SortLexicographic<2,T>());
	hull.clear();

	if(points.size()>=3){

	std::vector<Point<2,T>> & list = hull.getVertices();


	// Build lower hull
	for (int i = 0; i < points.size(); ++i) {
		while(list.size()>=2 && cross(list[list.size()-2],list[list.size()-1],points[i])<=0)
			list.pop_back();
		list.push_back(points[i]);
	}

	// Build upper hull
	for (int i =points.size()-2, t = list.size()+1; i >= 0; i--) {
		while (list.size() >= t && cross(list[list.size()-2],list[list.size()-1],points[i])<=0)
			list.pop_back();
		list.push_back(points[i]);
	}

	assert(list.size()==0 || ( list[0]==list.back()));
	if(list.size())
		list.pop_back();//for now, we aren't replicating the first vertex at the end of the polygon.
/*	for(auto & p:list){
		std::cout<<p.x << "," << p.y << " ";
	}
	std::cout<<"\n";*/
	hull.reorderVertices();
	}else{
		for(auto & p:points)
			hull.addVertex(p);
	}
	hull.update();
#ifndef NDEBUG
	for(int i = 0;i<pointSet.size();i++){
		Point<2,T> & p = pointSet[i];
		if(!pointSet.pointEnabled(i)){

		}else{
			assert(hull.contains(p));
		}
	}
#endif

}

template<>
inline mpq_class MonotoneConvexHull<2,mpq_class>::cross(const Point<2,mpq_class> &O, const Point<2,mpq_class> &A, const Point<2,mpq_class> &B)
{
	return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
}


#endif
