/*
 * MonotoneConvexHull.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */
#include "GeometryTypes.h"
#include "MonotoneConvexHull.h"
#include "mtl/Sort.h"
template<>
void MonotoneConvexHull<1,double>::update(){

}
template<>
void MonotoneConvexHull<2,double>::update(){
		vec<Point2D> points;
		pointSet.getEnabledPoints(points);
		sort(points,SortLexicographic<2,double>());
		hull.clear();
		vec<Point2D> & list = hull.getVertices();
		list.push(points[0]);
		list.push(points[1]);
		for (int i = 3;i<points.size();i++){

			while(list.size()>2){
				//check if the last three points do not make a right turn:
				if(cross(list[list.size()-2],list[list.size()-1],points[i])>=0){
					list.pop();
				}else{
					list.push(points[i]);
					break;
				}

			}
		}
		int cursize = list.size();
		for (int i = points.size()-2;i>=0;i--){
			while(list.size()>cursize){
				//check if the last three points do not make a right turn:
				if(cross(list[list.size()-2],list[list.size()-1],points[i])>=0){
					list.pop();
				}else{
					list.push(points[i]);
					break;
				}

			}
		}
		assert(list[0]==list.last());

	}

template<>
double MonotoneConvexHull<2,double>::getArea(){
	double area = 0;
	//traverse the polygons in clockwise order and compute the determinant
	//is there a better way to do this?
	int prevPID = -1;
	for (int i = 0;i<pointSet.getClockwisePoints().size();i++){
		int pid = pointSet.getClockwisePoints()[i];
		if(pointSet.pointEnabled(pid)){
			if(prevPID){
				auto & p = pointSet[pid];
				auto & prev = pointSet[prevPID];
				area+=prev.x *p.y - prev.y*p.x;
			}
			prevPID = pid;
		}
	}
	return area/2.0;
}

template<>
double MonotoneConvexHull<1,double>::getArea(){
	return 0;
}
