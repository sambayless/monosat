/*
 * MonotoneConvexHull.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */
#include "GeometryTypes.h"
#include "MonotoneConvexHull.h"
#include<algorithm>
#include "mtl/Sort.h"
template<>
void MonotoneConvexHull<1,double>::update(){

}
template<>
void MonotoneConvexHull<2,double>::update(){

		std::vector<Point2D> points;
		pointSet.getEnabledPoints(points);
		std::sort(points.begin(),points.end(),SortLexicographic<2,double>());
		hull.clear();
		std::vector<Point2D> & list = hull.getVertices();


		// Build lower hull
		for (int i = 0; i < points.size(); ++i) {
			if(list.size()>=2 && cross(list[list.size()-2],list[list.size()-1],points[i])<=0)
				list.pop_back();
			list.push_back(points[i]);
		}

		// Build upper hull
		for (int i =points.size()-2, t = list.size()+1; i >= 0; i--) {
			while (list.size() >= t && cross(list[list.size()-2],list[list.size()-1],points[i])<=0)
				list.pop_back();
			list.push_back(points[i]);
		}

		/*for (int i = 2;i<points.size();i++){

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
		}*/
		assert(list.size()==0 || ( list[0]==list.back()));
		if(list.size())
			list.pop_back();//for now, we aren't replicating the first vertex at the end of the polygon.
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
			if(prevPID>=0){
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
