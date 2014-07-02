/*
 * MonotoneDelaunay.h
 *
 *  Created on: Jun 30, 2014
 *      Author: sam
 */

#ifndef MONOTONE_DELAUNAY_H_
#define MONOTONE_DELAUNAY_H_
#include "Delaunay.h"

#include "PolygonSet.h"
#include <gmpxx.h>
#include <vector>

template<unsigned int D,class T>
class MonotoneDelaynay:public Delaunay<D,T>{

	PolygonSet<D,T> & polygons;

	struct MonotoneVertex{
		int polygonID;
		int pointID;
		vec<int> diagonals;//indices of diagonally connected vertices
		//int prev=-1;
		//int next=-1;
		bool seen=false;
		MonotoneVertex(int polygonID, int pointID,int nextID, int prevID):polygonID(polygonID),pointID(pointID),prev(prevID),next(nextID){};
	};
	enum VertexType{START,END,SPLIT,MERGE,REGULAR} ;
		struct Vertex{
			int polygonID;
			int pointID;
			int index;
			VertexType type;
		};
	std::vector<MonotoneVertex> monotonePolygons;
	std::vector<Vertex> sorted;
	std::vector<ConvexPolygon<D,T> > triangulation;
public:
	MonotoneDelaynay(PolygonSet<D,T> & p):polygons(p){

	}

	void update();

	std::vector<ConvexPolygon<D,T> > & getTriangulation(){
		update();
		return triangulation;
	}

private:
	void buildMonotonePolygons();
	void triangulate();
	void addDiagonal(int from, int to);
	inline bool convex(Point<D,T>& p1, Point<D,T>& p2, Point<D,T>& p3) {
		return ((p3.y-p1.y)*(p2.x-p1.x)-(p3.x-p1.x)*(p2.y-p1.y)) > 0;
	}


	inline bool IsInside(Point<D,T>& p1, Point<D,T>& p2, Point<D,T>& p3, Point<D,T> &p) {
		if(IsConvex(p1,p,p2)) return false;
		if(IsConvex(p2,p,p3)) return false;
		if(IsConvex(p3,p,p1)) return false;
		return true;
	}

	inline bool Below(Point<D,T> &p1, Point<D,T> &p2) {
		if(p1.y < p2.y) return true;
		else if(p1.y == p2.y) {
			if(p1.x < p2.x) return true;
		}
		return false;
	}

	T cross(const Point<D,T> &O, const Point<D,T> &A, const Point<D,T> &B);
};

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
//from http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
template<>
inline double MonotoneDelaynay<2,double>::cross(const Point2D &O, const Point2D &A, const Point2D &B)
{
	return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
}
template<>
inline mpq_class MonotoneDelaynay<2,mpq_class>::cross(const Point<2,mpq_class> &O, const Point<2,mpq_class> &A, const Point<2,mpq_class> &B)
{
	return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
}

template<>
void MonotoneDelaynay<2,double>::update();

template<>
void MonotoneDelaynay<2,mpq_class>::update();

template<>
void MonotoneDelaynay<2,mpq_class>::buildMonotonePolygons();

template<>
void MonotoneDelaynay<2,double>::buildMonotonePolygons();

#endif /* MONOTONEDELAUNAY_H_ */
