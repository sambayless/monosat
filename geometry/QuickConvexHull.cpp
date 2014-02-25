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
#include "mtl/Sort.h"

#include "cevans/quickhull3D.h"
#include "cevans/quickhull2D.h"
#include "cevans/zero.h"
template<>
 double zero<double>::val = 1E-10;

template<>
void QuickConvexHull<3,double>::update(){
	vec<Point3D> points;
	pointSet.getEnabledPoints(points);
	std::vector<Point3D> pointvec;
	for(int i =0;i<points.size();i++)
		pointvec.push_back(points[i]);
	quickhull3D<Point3D,double>  hull(pointvec);
	for(int i = 0;i<hull.boundary.size();i++){
		printf("%d ",hull.boundary[i]);
	}
	printf("\n");
}
template<>
void QuickConvexHull<2,double>::update(){
	vec<Point2D> points;
	pointSet.getEnabledPoints(points);
	std::vector<Point2D> pointvec;
	for(int i =0;i<points.size();i++)
		pointvec.push_back(points[i]);
	quickhull2D<Point2D,double>  hull(pointvec);
	for(int i = 0;i<hull.boundary.size();i++){
		printf("%d ",hull.boundary[i]);
	}
	printf("\n");
}

