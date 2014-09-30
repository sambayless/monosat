/****************************************************************************************[Solver.h]
The MIT License (MIT)

Copyright (c) 2014, Sam Bayless

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/


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
mpq_class cevans::zero<mpq_class>::val = mpq_class(0);
template<>
void QuickConvexHull<2,mpq_class>::update(){
	std::vector<Point<2,mpq_class>> points;
		pointSet.getEnabledPoints(points);
		if(points.size()<3){
			//edge case...
			hull.clear();
			for (auto & p:points)
				hull.addVertex(p);
		}else{
		std::sort(points.begin(),points.end(),SortLexicographic<2,mpq_class>());//not sure if this is required or not...

			 cevans::quickhull2D<Point<2,mpq_class>,mpq_class>  chull(points);
			hull.clear();
			for(int i = 0;i<chull.boundary.size();i++){
				//printf("%d ",chull.boundary[i]);
				hull.addVertex(points[chull.boundary[i]]);
			}
		}
	#ifndef NDEBUG
			for(int i = 0;i<pointSet.size();i++){
				Point<2,mpq_class> & p = pointSet[i];
				if(!pointSet.pointEnabled(i)){

				}else{
					assert(hull.contains(p,true));
				}
			}
	#endif

}



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
				assert(hull.contains(p,true));
			}
		}
#endif


}

