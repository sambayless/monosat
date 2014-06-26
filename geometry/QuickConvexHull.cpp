/*
 * QuickConvexHull.cpp
 *
 *  Created on: 2014-02-24
 *      Author: sam
 */


#include "QuickConvexHull.h"

#include <vector>
#include "GeometryTypes.h"
#include "MonotoneConvexHull.h"
#include <vector>
#include <algorithm>

#include "cevans/quickhull3D.h"
#include "cevans/quickhull2D.h"
#include "cevans/zero.h"
template<>
 double cevans::zero<double>::val = 1E-10;
template<>
void QuickConvexHull<1,double>::update(){

}
template<>
void QuickConvexHull<3,double>::update(){
	std::vector<Point3D> points;
	pointSet.getEnabledPoints(points);
	std::sort(points.begin(),points.end(),SortLexicographic<3,double>());
/*	std::vector<Point3D> pointvec;
	for(int i =0;i<points.size();i++)
		pointvec.push_back(points[i]);*/
	 cevans::quickhull3D<Point3D,double>  chull(points);
	for(int i = 0;i<chull.boundary.size();i++){
		printf("%d ",chull.boundary[i]);
	}
	printf("\n");

}
template<>
void QuickConvexHull<2,double>::update(){
	std::vector<Point2D> points;
	pointSet.getEnabledPoints(points);
	if(points.size()<3){
		//edge case...
		hull.clear();
		for (auto & p:points)
			hull.addVertex(p);
	}else{
	std::sort(points.begin(),points.end(),SortLexicographic<2,double>());//not sure if this is required or not...

		 cevans::quickhull2D<Point2D,double>  chull(points);
		hull.clear();
		for(int i = 0;i<chull.boundary.size();i++){
			//printf("%d ",chull.boundary[i]);
			hull.addVertex(points[chull.boundary[i]]);
		}
	}
#ifndef NDEBUG
		for(int i = 0;i<pointSet.size();i++){
			Point<2,double> & p = pointSet[i];
			if(!pointSet.pointEnabled(i)){

			}else{
				assert(hull.contains(p));
			}
		}
#endif


}
template<>
double QuickConvexHull<1,double>::getArea(){
	return 0;
}
template<>
double QuickConvexHull<2,double>::getArea(){
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

