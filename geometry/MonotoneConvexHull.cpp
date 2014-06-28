/*
 * MonotoneConvexHull.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */
#include "GeometryTypes.h"
#include "MonotoneConvexHull.h"
#include <algorithm>
#include "mtl/Sort.h"
#include <iostream>

template<>
void MonotoneConvexHull<1,double>::update(){

}
template<>
void MonotoneConvexHull<2,double>::update(){

		std::vector<Point2D> points;
		pointSet.getEnabledPoints(points);
		std::sort(points.begin(),points.end(),SortLexicographic<2,double>());
		hull.clear();

		if(points.size()>=3){

		std::vector<Point2D> & list = hull.getVertices();
/*		for(auto & p:points){
			std::cout << "("<<p.x << "," << p.y << "), ";
		}
		std::cout<<"\n";*/

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
		hull.reorderVertices();
/*		for(auto & p:list){
			std::cout << "("<<p.x << "," << p.y << "), ";
		}
		std::cout<<"\n";*/
		}else{
			for(auto & p:points)
				hull.addVertex(p);
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
void MonotoneConvexHull<2,mpq_class>::update(){
	std::vector<Point<2,mpq_class>> points;
		pointSet.getEnabledPoints(points);
		std::sort(points.begin(),points.end(),SortLexicographic<2,mpq_class>());
		hull.clear();

		if(points.size()>=3){

		std::vector<Point<2,mpq_class>> & list = hull.getVertices();


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
		for(auto & p:list){
			std::cout<<p.x << "," << p.y << " ";
		}
		std::cout<<"\n";
		hull.reorderVertices();
		}else{
			for(auto & p:points)
				hull.addVertex(p);
		}
#ifndef NDEBUG
		for(int i = 0;i<pointSet.size();i++){
			Point<2,mpq_class> & p = pointSet[i];
			if(!pointSet.pointEnabled(i)){

			}else{
				assert(hull.contains(p));
			}
		}
#endif
}


