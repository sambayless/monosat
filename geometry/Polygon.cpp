/*
 * Polygon.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#include "Polygon.h"
#include <cmath>
#include <algorithm>
//Note that this is subject to rounding errors. It also might only be correct for convex polygons.
//Also, this is only correct for planar (non-self-intersecting) polygons.
template<>
double Polygon<2,double>::getArea(){
	std::vector<Point2D> &  points = getVertices();

	double sum = 0;
	for (int i = 0;i<points.size();i++){
		Point2D& prev = i>0? points[i-1]: points.back();
		Point2D& cur = points[i];
		sum += prev[0]*cur[1]-cur[0]*prev[1];
	}
	return std::abs(sum/2.0);
}

template<>
bool Polygon<2,double>::contains(const Point<2,double> & point){
	std::vector<Point2D> &  points = getVertices();
	int i;
	  int j;
	  bool result = false;
	  for (i = 0, j = points.size() - 1; i < points.size(); j = i++) {
		if ((points[i].y > point.y) != (points[j].y > point.y) &&
			(point.x < (points[j].x - points[i].x) * (point.y - points[i].y) / (points[j].y-points[i].y) + points[i].x)) {
		  result = !result;
		 }
	  }
	  return result;
}

//Note that this is subject to rounding errors.
template<>
double Polygon<2,double>::getPerimeter(){
	std::vector<Point2D> & w = getVertices();
	double sum = 0;
	for (int i = 1;i<w.size();i++){
		Point2D prev = w[i-1];
		Point2D cur = w[i];
		double xdist = cur[0]-prev[0];
		double ydist=cur[1]-prev[1];
		sum += sqrt(xdist*xdist + ydist*ydist);
	}
	return sum;
}
template<>
void Polygon<2,double>::reorderVertices(){

	double centerX = 0;
	double centerY = 0;
	for(auto & p:vertices){
		centerX+=p.x;
		centerY+=p.y;
	}
	centerX/=vertices.size();
	centerY/=vertices.size();

	//from http://stackoverflow.com/a/6989383

	struct clockwise_lt{
		const std::vector<Point<2,double>> & points;
		double centerX;
		double centerY;
		clockwise_lt(const std::vector<Point<2,double>> & points,double centerX, double centerY):points(points),centerX(centerX),centerY(centerY){

		}

		bool operator () (int id_a, int id_b)
		{
			assert(id_a>=0);assert(id_a<points.size());
			assert(id_b>=0);assert(id_b<points.size());
			auto & a = points[id_a];
			auto & b = points[id_b];
			if (a[0] -centerX >= 0 && b[0] -centerX < 0)
				return true;
			if (a.x -centerX < 0 && b.x -centerX >= 0)
				return false;
			if (a.x -centerX == 0 && b.x -centerX == 0) {
				if (a.y - centerY >= 0 || b.y - centerY >= 0)
					return a.y > b.y;
				return b.y > a.y;
			}

			// compute the cross product of vectors (center -> a) x (center -> b)
			double det = (a.x -centerX) * (b.y - centerY) - (b.x -centerX) * (a.y - centerY);
			if (det < 0)
				return true;
			if (det > 0)
				return false;

			// points a and b are on the same line from the center
			// check which point is closer to the center
			double d1 = (a.x -centerX) * (a.x -centerX) + (a.y - centerY) * (a.y - centerY);
			double d2 = (b.x -centerX) * (b.x -centerX) + (b.y - centerY) * (b.y - centerY);
			return d1 > d2;
		}
	};
	//this should ideally be avoided...
	static std::vector<int> points_clockwise;
	points_clockwise.clear();
	for(int i =0;i<vertices.size();i++){
		points_clockwise.push_back(i);
	}
	std::sort(points_clockwise.begin(), points_clockwise.end(),clockwise_lt(vertices,centerX,centerY));
	//do this in place later
	static std::vector<Point<2,double>> oldPoints;
	oldPoints.clear();
	for(int i =0;i<vertices.size();i++){
		oldPoints.push_back(vertices[i]);
	}
	for(int i =0;i<vertices.size();i++){
		vertices[i] = oldPoints[points_clockwise[i]];
	}

	vertices_clockwise=true;
	assert(dbg_orderClockwise());
}
template<>
bool Polygon<1,double>::contains(const Point<1,double> & point){
	return false;
}
template<>
void Polygon<1,double>::reorderVertices(){

}
template<>
double Polygon<1,double>::getArea(){
	 return 0;
}
//Note that this is subject to rounding errors.
template<>
double Polygon<1,double>::getPerimeter(){
	return 0;
}

/*

//from http://content.gpwiki.org/index.php/Polygon_Collision
double determinant(Vector2D vec1, Vector2D vec2){
    return vec1.x * vec2.y - vec1.y * vec2.x;
}

//one edge is a-b, the other is c-d
Vector2D edgeIntersection(Vector2D a, Vector2D b, Vector2D c, Vector2D d){
    double det = determinant(b - a, c - d);
    double t   = determinant(c - a, c - d) / det;
    double u   = determinant(b - a, c - a) / det;
    if ((t < 0) || (u < 0) || (t > 1) || (u > 1)) {
        return NO_INTERSECTION;
    } else {
        return a * (1 - t) + t * b;
    }
}*/

template<>
void Polygon<2,mpq_class>::reorderVertices(){
	if(vertices.size()<3)
		return;
	mpq_class centerX = 0;
	mpq_class centerY = 0;
	for(auto & p:vertices){
		centerX+=p.x;
		centerY+=p.y;
	}
	centerX/=vertices.size();
	centerY/=vertices.size();

	//from http://stackoverflow.com/a/6989383

	struct clockwise_lt{
		const std::vector<Point<2,mpq_class>> & points;
		mpq_class centerX;
		mpq_class centerY;
		clockwise_lt(const std::vector<Point<2,mpq_class>> & points,mpq_class centerX, mpq_class centerY):points(points),centerX(centerX),centerY(centerY){

		}

		bool operator () (int id_a, int id_b)
		{
			assert(id_a>=0);assert(id_a<points.size());
			assert(id_b>=0);assert(id_b<points.size());
			auto & a = points[id_a];
			auto & b = points[id_b];
			if (a[0] -centerX >= 0 && b[0] -centerX < 0)
				return true;
			if (a.x -centerX < 0 && b.x -centerX >= 0)
				return false;
			if (a.x -centerX == 0 && b.x -centerX == 0) {
				if (a.y - centerY >= 0 || b.y - centerY >= 0)
					return a.y > b.y;
				return b.y > a.y;
			}

			// compute the cross product of vectors (center -> a) x (center -> b)
			auto det = (a.x -centerX) * (b.y - centerY) - (b.x -centerX) * (a.y - centerY);
			if (det < 0)
				return true;
			if (det > 0)
				return false;

			// points a and b are on the same line from the center
			// check which point is closer to the center
		 	auto d1 = (a.x -centerX) * (a.x -centerX) + (a.y - centerY) * (a.y - centerY);
			auto d2 = (b.x -centerX) * (b.x -centerX) + (b.y - centerY) * (b.y - centerY);
			return d1 > d2;
		}
	};
	//this should ideally be avoided...
	static std::vector<int> points_clockwise;
	points_clockwise.clear();
	for(int i =0;i<vertices.size();i++){
		points_clockwise.push_back(i);
	}
	std::sort(points_clockwise.begin(), points_clockwise.end(),clockwise_lt(vertices,centerX,centerY));
	//do this in place later
	static std::vector<Point<2,mpq_class>> oldPoints;
	oldPoints.clear();
	for(int i =0;i<vertices.size();i++){
		oldPoints.push_back(vertices[i]);
	}
	for(int i =0;i<vertices.size();i++){
		vertices[i] = oldPoints[points_clockwise[i]];
	}

	vertices_clockwise=true;
	assert(dbg_orderClockwise());
}

template<>
mpq_class Polygon<2,mpq_class>::getArea(){
	std::vector<Point<2,mpq_class>> &  points = getVertices();

	mpq_class sum = 0;
	for (int i = 0;i<points.size();i++){
		Point<2,mpq_class>& prev = i>0? points[i-1]: points.back();
		Point<2,mpq_class>& cur = points[i];
		sum += prev[0]*cur[1]-cur[0]*prev[1];
	}
	sum = sum/2;
	return abs(sum);
}

template<>
bool Polygon<2,mpq_class>::contains(const Point<2,mpq_class> & point){
	std::vector<Point<2,mpq_class>> &  points = getVertices();
	int i;
	  int j;
	  bool result = false;
	  for (i = 0, j = points.size() - 1; i < points.size(); j = i++) {
		if ((points[i].y > point.y) != (points[j].y > point.y) &&
			(point.x < (points[j].x - points[i].x) * (point.y - points[i].y) / (points[j].y-points[i].y) + points[i].x)) {
		  result = !result;
		 }
	  }
	  return result;
}

//Note that this is subject to rounding errors.
template<>
mpq_class Polygon<2,mpq_class>::getPerimeter(){
	std::vector<Point<2,mpq_class>> & w = getVertices();
	double sum = 0;
	/*for (int i = 1;i<w.size();i++){
		Point<2,mpq_class> prev = w[i-1];
		Point<2,mpq_class> cur = w[i];
		double xdist = cur[0]-prev[0];
		double ydist=cur[1]-prev[1];
		sum += sqrt(xdist*xdist + ydist*ydist);
	}*/
	return 0;
}
