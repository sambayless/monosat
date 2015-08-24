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

#ifndef CONVEXPOLYGON_H_
#define CONVEXPOLYGON_H_
#include "Polygon.h"
#include <gmpxx.h>
#include "core/Config.h"
#include "Line.h"
#include <cmath>
#include <iostream>
template<unsigned int D, class T>
class NConvexPolygon;

template<unsigned int D, class T>
class Triangle;

template<unsigned int D, class T>
class LineSegment;
/**
 * A concrete, convex polygon (or, for D>2, a polytope)
 */
template<unsigned int D, class T>
class ConvexPolygon: public Polygon<D, T> {
public:
	static long stats_triangle_avoided;
	static long stats_bounds_avoided;
	static long stats_split_checks;
	static long stats_contain_checks;
	static long stats_split_full_checks;
	static long stats_split_checks_depths;
	static long stats_bounds_intersections_avoided;
	static long stats_intersections;
	ConvexPolygon() :
			Polygon<D, T>() {
		
	}
	
	virtual ~ConvexPolygon() {
	}
	;
	virtual ShapeType getType() {
		return CONVEX_POLYGON;
	}
	bool contains(const Point<D, T> & point, bool inclusive);
	bool contains(const Point<D, T> & point, NConvexPolygon<D, T> * polygon_out, bool inclusive);

	bool containsInRange(const Point<D, T> & point, int firstVertex, int lastVertex, bool inclusive,
			bool excludeVertices = false);
	bool containsInSplit(const Point<D, T> & point, int firstVertex, int lastVertex, bool inclusive,
			bool excludeVertices = false);
	bool intersects(Shape<D, T> & s, bool inclusive);

	bool intersects(Shape<D, T> & s, NConvexPolygon<D, T> * polygon_out_this, NConvexPolygon<D, T> * polygon_out_other,
			bool inclusive);
	bool intersectsExcludingVertices(Shape<D, T> & s, NConvexPolygon<D, T> & polygon_out, bool inclusive);
	bool containsExcludingVertices(const Point<D, T> & point, NConvexPolygon<D, T> & polygon_out, bool inclusive);

	bool isConvex() {
		assert(dbg_Convex());
		return true;
	}
	
private:
	//bool intersects2d(Shape<2,T> & s, bool inclusive, bool ignore_vertices=false);
	bool intersects2d(Shape<2, T> & s, NConvexPolygon<2, T> * polygon_out_this,
			NConvexPolygon<2, T> * polygon_out_other, bool inclusive, bool ignore_vertices = false);
	bool edgesIntersectLine2d(LineSegment<2, T> & check, NConvexPolygon<2, T> * intersection, bool inclusive);
	bool containsInRange2d(const Point<2, T> & point, int firstVertex, int lastVertex, bool inclusive,
			bool excludeVertices);
	bool containsInSplit2d(const Point<2, T> & point, int firstVertex, int lastVertex,
			NConvexPolygon<2, T> & polygon_out, bool inclusive, bool excludeVertices = false);
	bool containsInSplit2d_helper(const Point<2, T> & point, int firstVertex, int lastVertex,
			NConvexPolygon<2, T> & polygon_out, int depth, bool inclusive, bool excludeVertices);
	bool dbg_Convex() {
#ifndef NDEBUG
		bool seenPositive = false;
		bool seenNegative = false;
		for (int i = 0; i < this->size(); i++) {
			const Point<2, T> & prev = i > 0 ? (*this)[i - 1] : (*this).back();
			const Point<2, T> & p = (*this)[i];
			const Point<2, T> & next = i < this->size() - 1 ? (*this)[i + 1] : (*this)[0];
			Point<2, T> a = p - prev;
			Point<2, T> b = next - p;
			T s = cross2d(a, b);
			seenPositive |= s > 0;
			seenNegative |= s < 0;
			if (seenPositive && seenNegative)
				return false;
		}
#endif
		return true;
		
	}
	
	static bool dbg_orderClockwise2dTri(Point<2, T> p1, Point<2, T> p2, Point<2, T> p3) {
#ifndef NDEBUG
		std::vector<Point<2, T>> points;
		points.push_back(p1);
		points.push_back(p2);
		points.push_back(p3);
		
		T sum = 0;
		for (int i = 0; i < points.size(); i++) {
			Point<2, T> & a = i > 0 ? points[i - 1] : points.back();
			Point<2, T> & b = points[i];
			sum += (b.x - a.x) * (b.y + a.y);
		}
		assert(sum >= 0);
		return sum >= 0;
		
#endif
		return true;
	}
};

template<unsigned int D, class T>
long ConvexPolygon<D, T>::stats_triangle_avoided = 0;
template<unsigned int D, class T>
long ConvexPolygon<D, T>::stats_bounds_avoided = 0;
template<unsigned int D, class T>
long ConvexPolygon<D, T>::stats_bounds_intersections_avoided = 0;
template<unsigned int D, class T>
long ConvexPolygon<D, T>::stats_split_checks = 0;
template<unsigned int D, class T>
long ConvexPolygon<D, T>::stats_contain_checks = 0;
template<unsigned int D, class T>
long ConvexPolygon<D, T>::stats_split_full_checks = 0;

template<unsigned int D, class T>
long ConvexPolygon<D, T>::stats_split_checks_depths = 0;

template<unsigned int D, class T>
long ConvexPolygon<D, T>::stats_intersections = 0;

/**
 * A triangle defined by three vertices
 */
template<unsigned int D, class T>
class Triangle: public ConvexPolygon<D, T> {
public:
	//A line segment is defined by two end-points that it passes through
	Point<D, T> a;
	Point<D, T> b;
	Point<D, T> c;

	Triangle() {
		
	}
	Triangle(const Point<D, T> & a, const Point<D, T> & b, const Point<D, T> & c) :
			a(a), b(b), c(c) {
		
	}
	int size() const {
		return 3;
	}
	
	void update() {
		this->bounds_uptodate = false;
	}
	
	const Point<D, T>& operator [](int index) const {
		index = index % size();
		if (index < 0) {
			index += size();
		}
		assert(index >= 0);
		assert(index < size());
		if (index == 0) {
			return a;
		} else if (index == 1) {
			return b;
		} else {
			assert(index == 2);
			return c;
		}
	}
	Point<D, T>& operator [](int index) {
		index = index % size();
		if (index < 0) {
			index += size();
		}
		assert(index >= 0);
		assert(index < size());
		if (index == 0) {
			return a;
		} else if (index == 1) {
			return b;
		} else {
			assert(index == 2);
			return c;
		}
	}
	
};

#include "LineSegment.h"

template<unsigned int D, class T>
class NConvexPolygon: public ConvexPolygon<D, T> {
public:
	
	std::vector<Point<D, T>> vertices;

	NConvexPolygon() :
			ConvexPolygon<D, T>() {
		
	}
	
	virtual ~NConvexPolygon() {
	}
	;

	int size() const {
		return vertices.size();
	}
	
	void update() {
		if (!this->vertices_clockwise) {
			reorderVertices();
		}
		this->bounds_uptodate = false;
	}
	
	void clear() {
		vertices.clear();
		this->bounds_uptodate = false;
	}
	
	void addVertex(Point<D, T> p) {
		this->vertices_clockwise = false;
		this->bounds_uptodate = false;
		vertices.push_back(p);
		this->distinct_vertices = -1;
	}
	/*	void addVertex(Point<D,T> & p){
	 vertices_clockwise=false;
	 bounds_uptodate=false;
	 vertices.push_back(p);
	 }*/
	//add a vertex, assuming that it will preserve clockwise order
	void addVertexUnchecked(Point<D, T> p) {
		if (vertices.size()) {
			if (p != vertices[0] && p != vertices.back()) {
				this->distinct_vertices++;
			}
		} else {
			this->distinct_vertices = 1;
		}
		vertices.push_back(p);
		assert(this->dbg_orderClockwise());
		assert(this->dbg_boundsUpToDate());
		assert(this->dbg_distinctVertexCount());
	}
	
	void popVertex() {
		vertices.pop_back();
		this->bounds_uptodate = false;
	}
	
	//Returns the vertices of the polygon, in clockwise order.
	
	std::vector<Point<D, T> > & getVertices() {
		if (!this->vertices_clockwise) {
			reorderVertices();
		}
		this->dbg_orderClockwise();
		return vertices;
	}
	const Point<D, T>& operator [](int index) const {
		index = index % size();
		if (index < 0) {
			index += size();
		}
		assert(index >= 0);
		assert(index < size());
		return vertices[index];
	}
	Point<D, T>& operator [](int index) {
		index = index % size();
		if (index < 0) {
			index += size();
		}
		assert(index >= 0);
		assert(index < size());
		return vertices[index];
	}
	
	//put the vertices into clockwise order
	void reorderVertices() {
		if (D == 2) {
			reorderVertices2d();
		} else
			assert(false);
	}
	
private:
	// copy ops are private to prevent copying
	//Polygon(const Polygon& from); // no implementation
	NConvexPolygon& operator=(const NConvexPolygon& from); // no implementation
			
	//put the vertices into clockwise order
	void reorderVertices2d();
	
};

template<unsigned int D, class T>
std::ostream & operator<<(std::ostream & str, ConvexPolygon<D, T> & polygon) {
	str << "ConvexPolygon=[";
	for (const auto & p : polygon) {
		str << p << ",";
	}
	str << "]";
	return str;
}

template<unsigned int D, class T>
std::ostream & operator<<(std::ostream & str, Triangle<D, T> & polygon) {
	str << "Triangle=[";
	for (const auto & p : polygon) {
		str << p << ",";
	}
	str << "]";
	return str;
}

template<unsigned int D, class T>
bool ConvexPolygon<D, T>::contains(const Point<D, T> & point, NConvexPolygon<D, T> * polygon_out, bool inclusive) {
	if (polygon_out)
		polygon_out->clear();
	stats_contain_checks++;
	if (this->size() == 3) {
		if (D == 2) {
			stats_triangle_avoided++;
			if (this->pointInTriangle2d(point, (*this)[0], (*this)[1], (*this)[2], inclusive)) {
				if (polygon_out) {
					polygon_out->addVertexUnchecked((*this)[0]);
					polygon_out->addVertexUnchecked((*this)[1]);
					polygon_out->addVertexUnchecked((*this)[2]);
				}
				return true;
			} else {
				return false;
			}
		}
	}
	if (!this->boundContains(point, inclusive)) {
		assert(!containsInRange(point, 0, this->size() - 1, inclusive));
		stats_bounds_avoided++;
		return false;
	}
	
	if (D == 2) {
		if (polygon_out)
			return containsInSplit2d((const Point<2, T> &) point, 0, this->size() - 1, *polygon_out, inclusive);
		else {
			NConvexPolygon<2, T> ignore;
			return containsInSplit2d((const Point<2, T> &) point, 0, this->size() - 1, ignore, inclusive);
		}
	} else {
		assert(false);
	}
	return false;
}

template<unsigned int D, class T>
bool ConvexPolygon<D, T>::contains(const Point<D, T> & point, bool inclusive) {
	return contains(point, nullptr, inclusive);
	/*	//stats_contain_checks++;
	 static int iter = 0;
	 if(++iter==309222){
	 int a=1;
	 }
	 if(this->size()==1){
	 if(inclusive){
	 return (*this)[0]==point;
	 }else
	 return false;
	 }
	 if( this->size()==3){
	 if(D== 2)
	 return this->pointInTriangle2d(point,(*this)[0],(*this)[1],(*this)[2],inclusive);
	 }
	 if(!this->boundContains(point,inclusive)){
	 //bool b = Polygon<D,T>::boundContains(point);
	 assert(!containsInRange(point,0,this->size()-1,inclusive));
	 //stats_bounds_avoided++;
	 return false;
	 }

	 if(Minisat::pipalg == PointInPolygonAlg::ALG_RECURSIVE_SPLIT){

	 return containsInSplit(point,0,this->size()-1,inclusive);
	 }else{
	 return containsInRange(point,0,this->size()-1,inclusive);
	 }*/
}
template<unsigned int D, class T>
bool ConvexPolygon<D, T>::containsInRange(const Point<D, T> & point, int firstVertex, int lastVertex, bool inclusive,
		bool excludeVertices) {
	this->update();
	if (D == 2) {
		return containsInRange2d((const Point<2, T> &) point, firstVertex, lastVertex, inclusive, excludeVertices);
	} else {
		assert(false);
	}
	return false;
}
template<unsigned int D, class T>
bool ConvexPolygon<D, T>::containsInSplit(const Point<D, T> & point, int firstVertex, int lastVertex, bool inclusive,
		bool excludeVertices) {
	this->update();
	if (D == 2) {
		static NConvexPolygon<2, T> ignore;
		return containsInSplit2d((const Point<2, T> &) point, firstVertex, lastVertex, ignore, inclusive,
				excludeVertices);
	} else {
		assert(false);
	}
	return false;
}

template<unsigned int D, class T>
bool ConvexPolygon<D, T>::containsInSplit2d(const Point<2, T> & point, int firstVertex, int lastVertex,
		NConvexPolygon<2, T> & polygon_out, bool inclusive, bool excludeVertices) {
	stats_split_checks++;
	//std::vector<Point<2,T> > &  w = (std::vector<Point<2,T> > & ) Polygon<D,T>::getVertices();
	ConvexPolygon<2, T> & w = (ConvexPolygon<2, T> &) *this;
	polygon_out.clear();
	if (w.size() == 0)
		return false;
	if (lastVertex < 0)
		lastVertex = w.size() - 1;
	assert(lastVertex >= 0);
	int n_verts = lastVertex - firstVertex + 1;
	if (lastVertex < firstVertex) {
		n_verts = lastVertex + (w.size() - firstVertex);
	}
	assert(n_verts <= w.size());
	this->dbg_orderClockwise();
	if (n_verts == 0)
		return false;
	else if (n_verts == 1) {
		if (!inclusive || excludeVertices)
			return false;
		assert(lastVertex == firstVertex);
		if (w[firstVertex] == point) {
			polygon_out.addVertex(w[firstVertex]);
			return true;
		}
		return false;
	} else if (n_verts == 2) {
		if (!inclusive)
			return false;
		assert(lastVertex != firstVertex);
		//from http://stackoverflow.com/a/11908158
		//true if the point is between (inclusive) the other two points.
		auto p1 = w[firstVertex];
		auto p2 = w[lastVertex];
		
		if (excludeVertices) {
			if (p1 == point || p2 == point) {
				return false;
			}
		}
		
		//check if the point lies on this line
		if (crossDif(point, p1, p2) == 0) {
			T dxl = p2.x - p1.x;
			T dyl = p2.y - p1.y;
			bool contains;
			using namespace std;
			//check if the point is between the end points
			if (abs(dxl) >= abs(dyl))
				contains = dxl > 0 ? p1.x <= point.x && point.x <= p2.x : p2.x <= point.x && point.x <= p1.x;
			else
				contains = dyl > 0 ? p1.y <= point.y && point.y <= p2.y : p2.y <= point.y && point.y <= p1.y;
			
			if (contains) {
				polygon_out.addVertex(w[firstVertex]);
				polygon_out.addVertex(w[lastVertex]);
			}
			return contains;
		}
		return false;
	} else if (n_verts == 3) {
		bool contains = containsInRange(point, firstVertex, lastVertex, inclusive, excludeVertices);
		if (contains) {
			polygon_out.addVertex(w[firstVertex]);
			polygon_out.addVertex(w[firstVertex + 1]);
			polygon_out.addVertex(w[firstVertex + 2]);
		}
		return contains;
	}
	stats_split_full_checks++;
	
	bool res = containsInSplit2d_helper(point, firstVertex, lastVertex, polygon_out, 1, inclusive, excludeVertices);
	if (res && !inclusive) {
		//do an extra check to correct for edge cases.
		assert(polygon_out.contains(point, true));
		if (!polygon_out.contains(point, false)) {
			//then there are two possibilities. Either the point is exactly on an edge of the polygon (in whch case it is not contained unless we are being inclusive)
			//OR, we can expand the triangle into a quadrilateral that _does_ contain the point.
			//so, check which of the three edges of the triangle the point is lying on (if it is on a vertex of the triangle, then it must be of the border of the polygon, so it is not included)
			assert(polygon_out.size() == 3);
			int containing_edge = -1;
			for (int i = 0; i < polygon_out.size(); i++) {
				Point<2, T> & prev = polygon_out[i - 1];
				Point<2, T> & p = polygon_out[i];
				if (eq_epsilon(crossDif(prev, p, point))) {
					containing_edge = i;
					//check if it is on a vertex exactly
					if (point == p || point == prev) {
						//point is exactly on a vertex, so is not included
						polygon_out.clear();
						return false;
					}
					
					break;
				}
			}
			assert(containing_edge >= 0);
			int pIndex = -1;
			int prevIndex = -1;
			int nIndex = -1;
			//now check whether that edge is on the border of the polygon
			for (int i = 0; i < this->size(); i++) {
				Point<2, T> & prev = (*this)[i - 1];
				Point<2, T> & p = (*this)[i];
				if (p == polygon_out[containing_edge]) {
					pIndex = i;
				} else if (p == polygon_out[containing_edge - 1]) {
					prevIndex = i;
				} else if (p == polygon_out[containing_edge + 1]) {
					nIndex = i;
				}
			}
			if (pIndex < prevIndex) {
				std::swap(pIndex, prevIndex);
			}
			
			if (((pIndex + 1) % this->size() == prevIndex || (pIndex - 1) % this->size() == prevIndex)
					|| ((*this)[pIndex + 1] == (*this)[prevIndex] || (*this)[pIndex - 1] == (*this)[prevIndex])) {
				//The containing edge is on the border of the polygon, so the point is NOT contained exclusively within the polygon.
				//Note that there is a corner case here, if there are two or more IDENTICAL points in the convex hull. In this case, it is possible for the
				//point to be exactly on an edge, but for the indices of those points (in clockwise order) to not be adjacent.
				//So it really is necessary to check the actual values of the points in the above check.
				res = false;
				polygon_out.clear();
				
			} else {
				
				if (res) {
					//None of the edges of the hull contain the point, so the containing edge is interior to the polygon - meaning that the point is exclusively contained in the polygon, but that the current containing triangle
					//either needs to be swapped for some other triangle, or widened to a quadrilateral.
					//(since there exist edge cases where NO internal triangle will exclusively contain the point, we will just widen to a quadrilateral, here.)
					assert(prevIndex < pIndex);
					assert(pIndex > -1);
					assert(prevIndex > -1);
					assert(nIndex > -1);
					assert(pIndex != prevIndex);
					assert(pIndex != nIndex);
					
					if (prevIndex < nIndex && nIndex < pIndex) {
						//ok, pick any vertex outside of pIndex, prevIndex in order to expand the quadrilateral to contain point
						assert((nIndex + 1) % this->size() != prevIndex);
						//We *should* be able to pick any vertex outside that range -
						//However, there is a corner case when the convex hull contains two identical points, in which case we don't get containment.
						int newIndex = pIndex + 1;
						
						while ((*this)[newIndex] == (*this)[pIndex]) {
							newIndex = (newIndex + 1) % this->size();
						}
						assert((*this)[newIndex] != (*this)[pIndex]);
						assert((*this)[newIndex] != (*this)[prevIndex]);
						polygon_out.addVertex((*this)[newIndex]);
					} else {
						assert(nIndex < prevIndex || nIndex > pIndex);
						assert((pIndex - 1) % this->size() != prevIndex);
						
						int newIndex = pIndex - 1;
						
						while ((*this)[newIndex] == (*this)[pIndex]) {
							newIndex = (newIndex - 1) % this->size();
						}
						assert((*this)[newIndex] != (*this)[pIndex]);
						assert((*this)[newIndex] != (*this)[prevIndex]);
						polygon_out.addVertex((*this)[newIndex]);
					}
					/*	std::cout<<*this<<"\n";
					 std::cout << polygon_out<< "\n";
					 std::cout <<point<<"\n";*/

					assert(polygon_out.containsInRange(point, 0, polygon_out.size() - 1, true));
					assert(polygon_out.containsInRange(point, 0, polygon_out.size() - 1, false));
				}
			}
			
		}
	}
	
	assert(res == containsInRange(point, firstVertex, lastVertex, inclusive));
	return res;
	
}
template<unsigned int D, class T>
bool ConvexPolygon<D, T>::containsInSplit2d_helper(const Point<2, T> & point, int first_vertex, int last_vertex,
		NConvexPolygon<2, T> & triangle_out, int depth, bool inclusive, bool excludeVertices) {
	//recurse on this segment of the polygon, finding a triangle that contains the point.
	//precondition: the point is contained in this convex segment of the polygon
	
	//Noah's algorithm: pick 3 vertices in the polygon. 2 of them are adjacent, and the third is arbitrary (but should probably be the index that is farthest in both directions from the adjacent vertices)
	//Check if they contain the point; if they do, return them.
	//Else, check which of the two sides of the triangle with the non-adjacent vertex the point is on. Recurse on that sub polygon.
	
	//When recursing, 2 of the three vertices are already selected (they are the vertices from the existing triangle), so we only have to pick one more vertex.
	//Since we already know that the point isn't on the other side of those two vertices, we only have to check two sides in the case where the point is not contained.
	assert(first_vertex != last_vertex);
	
	//std::vector<Point<2,T>> & polygon_vertices =(std::vector<Point<2,T>> & ) Polygon<D,T>::vertices;
	ConvexPolygon<2, T> & polygon_vertices = (ConvexPolygon<2, T> &) *this;
	Point<2, T> & a = polygon_vertices[first_vertex];
	Point<2, T> & b = polygon_vertices[last_vertex];
	assert(first_vertex < last_vertex);
	if (first_vertex + 1 == last_vertex) {
		stats_split_checks_depths += depth;
		if (!inclusive || excludeVertices) {
			return false;
		}
		//Then this is a line. depending on our notion of containment, we either give up, or test if the line contains this point
		assert(!excludeVertices);
		/*	if (excludeVertices && (a == point || b == point)) {
			//unreachable, because of above.
			return false;
		}*/

		triangle_out.clear();
		triangle_out.addVertexUnchecked(a);
		triangle_out.addVertexUnchecked(b);
		if (triangle_out.contains(point, inclusive)) {
			return true;
		} else {
			triangle_out.clear();	 //give up
		}
		assert(!containsInRange(point, first_vertex, last_vertex, inclusive));
		return false;
	}
	
	int mid_point = 0;
	if (first_vertex < last_vertex) {
		mid_point = (last_vertex - first_vertex) / 2 + first_vertex;
	} else {
		mid_point = (first_vertex - last_vertex) / 2 + last_vertex;
	}
	Point<2, T> & c = polygon_vertices[mid_point];
	assert(mid_point != last_vertex);
	assert(mid_point != first_vertex);
	/**/

	//compute the barycentric coordinates of this point
	auto & p0 = c;	 //b;
	auto & p1 = a;	 //c;
	auto & p2 = b;	 //a;
	
	assert(dbg_orderClockwise2dTri(p2, p1, p0));//intentionally reversing winding here, because the formula below is counter clockwise.
	T s = (p0.y * p2.x - p0.x * p2.y + (p2.y - p0.y) * point.x + (p0.x - p2.x) * point.y);
	T t = (p0.x * p1.y - p0.y * p1.x + (p0.y - p1.y) * point.x + (p1.x - p0.x) * point.y);
	T area2 = (-p1.y * p2.x + p0.y * (-p1.x + p2.x) + p0.x * (p1.y - p2.y) + p1.x * p2.y);
	assert(area2 >= 0);
	/*T s = (a.y*b.x - a.x*b.y + (b.y - a.y)*point.x + (a.x - b.x)*point.y);
	 T t = (a.x*c.y - a.y*c.x + (a.y - c.y)*point.x + (c.x - a.x)*point.y);*/
	//T gamma = 1-s-t;
	bool contained;
	
	//checks for inclusiveness will be done outside this method
	//if(inclusive)
	contained = (s >= 0 && t >= 0 && (s + t <= area2));
	//else
	//	contained= ( s>0 && t>0 && (s+t<area2));
	
	if (contained) {
		triangle_out.clear();
		
		if (excludeVertices && (a == point || b == point || c == point)) {
			return false;
		}
		
		triangle_out.addVertexUnchecked(a);
		triangle_out.addVertexUnchecked(c);
		triangle_out.addVertexUnchecked(b);
		assert(triangle_out.contains(point, true));
		stats_split_checks_depths += depth;
		return true;	//we are done
	} else {
#ifndef NDEBUG
		NConvexPolygon<D, T> dbg_poly;
		dbg_poly.addVertexUnchecked(a);
		dbg_poly.addVertexUnchecked(c);
		dbg_poly.addVertexUnchecked(b);
		assert(!dbg_poly.contains(point, true));
		
#endif
		
		if (t < 0 && s >= 0) {
			//point is either between first_vertex,mid_point, or not in the triangle
			return containsInSplit2d_helper(point, first_vertex, mid_point, triangle_out, depth + 1, inclusive,
					excludeVertices);
		} else if (s < 0 && t >= 0) {
			//point is either between first_vertex,mid_point, or not in the triangle
			return containsInSplit2d_helper(point, mid_point, last_vertex, triangle_out, depth + 1, inclusive,
					excludeVertices);
		} else {
			stats_split_checks_depths += depth;
			//point is not contained.
			assert(!containsInRange(point, first_vertex, last_vertex, inclusive));
			return false;	//this point is not contained.
		}
	}
}

template<unsigned int D, class T>
bool ConvexPolygon<D, T>::containsInRange2d(const Point<2, T> & point, int firstVertex, int lastVertex, bool inclusive,
		bool excludeVertices) {
	//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
	//this is correct _only_ for convex polygons
	//std::vector<Point<2,T> > &  w = (std::vector<Point<2,T> > & ) Polygon<D,T>::getVertices();
	ConvexPolygon<2, T> & w = (ConvexPolygon<2, T> &) *this;
	if (w.size() == 0)
		return false;
	if (lastVertex < 0)
		lastVertex = w.size() - 1;
	assert(lastVertex >= 0);
	int n_verts = lastVertex - firstVertex + 1;
	if (lastVertex < firstVertex) {
		n_verts = lastVertex + (w.size() - firstVertex);
	}
	assert(n_verts <= w.size());
	this->dbg_orderClockwise();
	if (n_verts == 0)
		return false;
	else if (n_verts == 1) {
		assert(lastVertex == firstVertex);
		return inclusive && (!excludeVertices) && w[firstVertex] == point;
	} else if (n_verts == 2) {
		if (!inclusive)
			return false;
		assert(lastVertex != firstVertex);
		//from http://stackoverflow.com/a/11908158
		//true if the point is between (inclusive) the other two points.
		auto p1 = w[firstVertex];
		auto p2 = w[lastVertex];
		
		if (excludeVertices) {
			if (p1 == point || p2 == point) {
				return false;
			}
		}
		//check if the point lies on this line
		if (crossDif(point, p1, p2) == 0) {
			T dxl = p2.x - p1.x;
			T dyl = p2.y - p1.y;
			bool contains;
			using namespace std;
			//check if the point is between the end points
			if (abs(dxl) >= abs(dyl)) {
				contains = dxl > 0 ? p1.x <= point.x && point.x <= p2.x : p2.x <= point.x && point.x <= p1.x;
			} else {
				contains = dyl > 0 ? p1.y <= point.y && point.y <= p2.y : p2.y <= point.y && point.y <= p1.y;
			}
			
			return contains;
		}
		return false;
	}
	
	int endVertex = (lastVertex + 1) % w.size();
	//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
	//this is correct _only_ for convex polygons
	//note: this can also compute the area (which is the sum of p2[0]*p1[1] - p1[0]*p2[1]); could potentially combine these...
	
	if (excludeVertices) {
		//check if the point is exactly a vertex, in which case the point is not considered to be contained.
		//can we do this in log time?
		for (int n = 0; n < n_verts; n++) {
			int i = (firstVertex + n) % w.size();
			if (w[i] == point) {
				return false;
			}
		}
	}
	
	for (int n = 0; n < n_verts; n++) {
		int i = (firstVertex + n) % w.size();
		Point<2, T> p1 = (i > 0 ? w[i - 1] : w.back()) - point;
		Point<2, T> p2 = w[i] - point;
		
		bool contained;
		if (inclusive)
			contained = (p2[0] * p1[1] - p1[0] * p2[1]) >= 0;
		else
			contained = (p2[0] * p1[1] - p1[0] * p2[1]) > 0;
		if (!contained) {
			return false;
		}
	}
	return true;
}

template<unsigned int D, class T>
bool ConvexPolygon<D, T>::containsExcludingVertices(const Point<D, T> & point, NConvexPolygon<D, T> & polygon_out,
		bool inclusive) {
	polygon_out.clear();
	stats_contain_checks++;
	
	if (!this->boundContains(point, inclusive)) {
		assert(!containsInRange(point, 0, this->size() - 1, inclusive));
		stats_bounds_avoided++;
		return false;
	}
	
	if (D == 2) {
		return containsInSplit2d((const Point<2, T> &) point, 0, this->size() - 1, (NConvexPolygon<D, T> &) polygon_out,
				inclusive, true);
	} else {
		assert(false);
	}
	return false;
}

template<unsigned int D, class T>
bool ConvexPolygon<D, T>::intersectsExcludingVertices(Shape<D, T> & shape, NConvexPolygon<D, T> & out, bool inclusive) {
	out.clear();
	if (D == 2) {
		return intersects2d((Shape<2, T> &) shape, (NConvexPolygon<2, T> &) out, inclusive, true);
	} else {
		assert(false);
	}
	return false;
}

template<unsigned int D, class T>
bool ConvexPolygon<D, T>::intersects(Shape<D, T> & shape, NConvexPolygon<D, T> * polygon_out_this,
		NConvexPolygon<D, T> * polygon_out_other, bool inclusive) {
	if (polygon_out_this) {
		polygon_out_this->clear();
	}
	if (polygon_out_other) {
		polygon_out_other->clear();
	}
	if (D == 2) {
		return intersects2d((Shape<2, T> &) shape, polygon_out_this, polygon_out_other, inclusive);
	} else {
		assert(false);
	}
	return false;
}
template<unsigned int D, class T>
bool ConvexPolygon<D, T>::intersects(Shape<D, T> & shape, bool inclusive) {
	if (D == 2) {
		return intersects2d((Shape<2, T> &) shape, nullptr, nullptr, inclusive);
	} else {
		assert(false);
	}
	return false;
}

template<unsigned int D, class T>
bool ConvexPolygon<D, T>::edgesIntersectLine2d(LineSegment<2, T> & check, NConvexPolygon<2, T> * out, bool inclusive) {
	// std::vector<Point<2,T> > &  w = this->getVertices();
	ConvexPolygon<2, T> & w = (ConvexPolygon<2, T>&) *this;
	assert(check.a != check.b);
	if (!inclusive && w.size() > 2) {
		//There is an edge case here where both end points of the line are exactly on vertices or edges of the polygon,
		//and the line also crosses through the polygon.
		//to check whether this occurs, pick any point on the line that is not an endpoint, and see if they collide.
		Point<2, T> mid = (check.a + check.b) / 2;
		if (this->contains(mid, out, inclusive)) {
			return true;
		}
	}
	if (out)
		out->clear();
	static LineSegment<2, T> test;
	bool seenContainedLine = false;
	for (int i = 0; i < w.size(); i++) {
		Point<2, T> & prev = i > 0 ? w[i - 1] : w.back();
		Point<2, T> & p = w[i];
		test.a = prev;
		test.b = p;
		//Have to be VERY careful about inclusive, here. A point that crosses through the polygon by passing through a vertex may still intersect the shape even if it collision detection is not inclusive!
		if (check.intersects(test, inclusive)) {
			if (out) {
				out->addVertex(prev);
				out->addVertex(p);
			}
			return true;
		} else if (!inclusive) {
			//we need to check if the endpoint is contained inside the line (exclusively of the endpoints)
			
			/*	if(p != check.a && p != check.b && check.contains(p,true)){
			 //this is a non-inclusive collision UNLESS this polygon has <= two distinct points.
			 if(this->distinctVertices()>2){

			 if(out){
			 out->addVertex(p);
			 }
			 return true;
			 }
			 }*/
			//there is also the case where both ends of check are vertices of the polygon, but the line still passes through the interior of the polygon.
			//that is handled elsewhere in the collision code.
		}
	}
	return false;
}
/*

 template<unsigned int D,class T>
 bool ConvexPolygon<D,T>::intersects2d(Shape<2,T> & shape, bool inclusive, bool ignore_vertices){
 if(this->size()==0)
 return false;
 if(shape.getType()==LINE_SEGMENT){
 LineSegment<2,T> & line = (LineSegment<2,T> &)shape;

 //first, check if either end point is contained
 if(this->contains(line.a,inclusive) || this->contains(line.b,inclusive))
 return true;
 if(this->size()==1){
 if(ignore_vertices)
 return false;
 return line.contains((*this)[0],inclusive);
 }
 static NConvexPolygon<2,T>  ignore;
 //the line may still intersect even if neither end point is contained.
 //we could apply the SAT here. But instead, we're going to walk around the edges of the convex shape, and see if any of the edges intersect this line.
 return edgesIntersectLine2d(line,ignore, inclusive);


 }else if(shape.getType()==CONVEX_POLYGON){
 ConvexPolygon<2,T> & c = (ConvexPolygon<2,T>&) shape;
 if(c.size()<this->size()){
 return c.intersects(*this,inclusive);
 }

 if(this->size()==0 || c.size()==0){
 return false;
 }else if (this->size()==1){
 return c.contains((*this)[0],inclusive);
 if(!inclusive){
 return false;
 }
 if(c.size()==1){
 //then the two vertices are considered to collide if they are identical
 return (*this)[0]==c[0];
 }else{
 return c.contains((*this)[0],inclusive);
 }
 }else if (c.size()==1){
 if(!inclusive)
 return false;
 return contains(c[0],inclusive);
 }


 ConvexPolygon<2,T> &  w = (ConvexPolygon<2,T>&)*this;

 //Separating Axis Theorem for collision detection between two convex polygons
 //loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
 //If any of the projections are non-intersection, then these don't collide; else, they do collide
 if(this->size()>1){
 for(int i = 0;i<this->size();i++){
 auto & cur = (*this)[i];
 auto & prev = (*this)[i-1];
 Point<2,T> edge = cur-prev;
 Point<2,T> un_normalized_normal(-edge.y, edge.x);

 //now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
 //Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
 //(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
 T left = numeric<T>::infinity();
 T right = -numeric<T>::infinity();
 for (auto & p:*this){
 T projection = un_normalized_normal.dot(p);
 if (projection < left) {
 left = projection;
 }
 if (projection > right) {
 right = projection;
 }
 }

 bool seenLeft = false;
 bool seenRight=false;
 for (auto & p:c){
 T projection = un_normalized_normal.dot(p);
 if(inclusive){
 if (projection >= left && projection <= right ) {
 seenRight=true;
 seenLeft=true;
 break;
 }else if (projection < left ){
 seenLeft=true;
 if(seenRight){
 break;
 }
 }else if (projection>right){
 seenRight=true;
 if (seenLeft){
 break;
 }
 }else if (seenLeft && projection > left ){
 seenRight=true;
 break;
 }else if (seenRight && projection < right ){
 seenRight=true;
 break;
 }
 }else{
 if(projection>left){
 seenRight=true;
 if (seenLeft){
 break;
 }
 }
 if (projection<right){
 seenLeft=true;
 if(seenRight){
 break;
 }
 }

 }
 }
 if(!(seenLeft&&seenRight)){
 return false;
 }
 }
 }
 if(c.size()>1){
 //now test the axis produced by the other polygon
 for(int j = 0;j<c.size();j++){
 auto & cur = c[j];
 auto & prev =c[j-1];
 Point<2,T> edge =cur-prev;
 Point<2,T> un_normalized_normal(-edge.y, edge.x);

 T left = numeric<T>::infinity();
 T right = -numeric<T>::infinity();
 for (auto & p:c){
 T projection = un_normalized_normal.dot(p);
 if (projection < left) {
 left = projection;
 }
 if (projection > right) {
 right = projection;
 }
 }
 bool seenLeft = false;
 bool seenRight=false;

 for (auto & p:*this){
 T projection = un_normalized_normal.dot(p);
 if(inclusive){
 if (projection >= left && projection <= right ) {
 seenRight=true;
 seenLeft=true;
 break;
 }else if (projection < left ){
 seenLeft=true;
 if(seenRight){
 break;
 }
 }else if (projection>right){
 seenRight=true;
 if (seenLeft){
 break;
 }
 }else if (seenLeft && projection > left ){
 seenRight=true;
 break;
 }else if (seenRight && projection < right ){
 seenRight=true;
 break;
 }
 }else{
 if(projection>left){
 seenRight=true;
 if (seenLeft){
 break;
 }
 }
 if (projection<right){
 seenLeft=true;
 if(seenRight){
 break;
 }
 }

 }
 }
 if(!(seenLeft&&seenRight)){
 return false;
 }
 }
 //If no axis overlapped, then they did in fact intersect
 return true;
 }
 }
 assert(false);
 return false;
 }
 */

template<unsigned int D, class T>
bool ConvexPolygon<D, T>::intersects2d(Shape<2, T> & shape, NConvexPolygon<2, T> * polygon_out_this,
		NConvexPolygon<2, T> * polygon_out_other, bool inclusive, bool ignore_vertices) {
	static int iter = 0;
	if (++iter == 3) {
		int a = 1;
	}
	if (this->size() == 0)
		return false;
	
	if (!this->boundIntersects(shape, inclusive)) {
#ifndef NDEBUG
		NConvexPolygon<2, T> t;
		for (auto & p : *this) {
			t.addVertex(p);
		}
		assert(!t.intersects2d(shape, nullptr, nullptr, inclusive, ignore_vertices));
#endif
		stats_bounds_intersections_avoided++;
		return false;
	}
	stats_intersections++;
	if (shape.getType() == LINE_SEGMENT) {
		LineSegment<2, T> & line = (LineSegment<2, T> &) shape;
		
		if (line.a == line.b) {
			//treat this line segment as a point.
			if (this->contains(line.a, polygon_out_this, inclusive)) {
				if (polygon_out_other) {
					polygon_out_other->clear();
					polygon_out_other->addVertex(line[0]);
				}
				return true;
			}
			return false;
		}
		
		//first, check if either end point is contained
		if (this->contains(line.a, polygon_out_this, inclusive)
				|| this->contains(line.b, polygon_out_this, inclusive)) {
			if (polygon_out_other) {
				polygon_out_other->clear();
				polygon_out_other->addVertex(line[0]);
				polygon_out_other->addVertex(line[1]);
			}
			return true;
		}
		if (this->size() == 1) {
			if (ignore_vertices)
				return false;
			if (line.contains((*this)[0], inclusive)) {
				if (polygon_out_this) {
					polygon_out_this->clear();
					polygon_out_this->addVertex((*this)[0]);
				}
				if (polygon_out_other) {
					polygon_out_other->clear();
					polygon_out_other->addVertex(line[0]);
					polygon_out_other->addVertex(line[1]);
				}
				return true;
			}
			return false;
		}
		
		//the line may still intersect even if neither end point is contained.
		//we could apply the SAT here. But instead, we're going to walk around the edges of the convex shape, and see if any of the edges intersect this line.
		//Have to be VERY careful about inclusive, here. A point that crosses through the polygon by passing through a vertex may still intersect the shape even if it is not inclusive!
		
		bool r = edgesIntersectLine2d(line, polygon_out_this, inclusive);
		if (r && polygon_out_other) {
			polygon_out_other->clear();
			polygon_out_other->addVertex(line[0]);
			polygon_out_other->addVertex(line[1]);
		}
		
		if (!inclusive && !r && this->size() >= 3) {
			//If collisions are not inclusive, AND no edge was intersected, it is still possible that the line directly passes through (or ends at)
			//two vertices of the polygon, or two edges of the polygon, so long as this polygon isn't a triangle (note that it can still be a triangle if it has more than 3 points, but some of them are exactly on an edge or are identical).
			if (polygon_out_this) {
				polygon_out_this->clear();
			}
			
			//First, if the line's endpoints lie on two different edges of the polygon, then the midpoint of the line must be contained.
			Point<2, T> midpoint = (line.a + line.b) / 2;
			if (this->contains(midpoint, polygon_out_this, inclusive)) {
				if (polygon_out_other) {
					polygon_out_other->clear();
					polygon_out_other->addVertex(line[0]);
					polygon_out_other->addVertex(line[1]);
				}
				return true;
			}
			
			//ok, then the line must pass through atleast one vertex of the polygon. The other end may ALSO pass through a vertex, or may land directly on an edge.
			
			//so first, find a vertex of the polygon that is on the line.
			int first_pass_vertex = -1;
			
			for (int i = 0; i < this->size(); i++) {
				Point<2, T> & p = (*this)[i];
				if (line.contains(p, true)) {			//this check is intentionally inclusive
					if (first_pass_vertex == -1) {
						first_pass_vertex = i;
					} else {
						//If the line passes through two vertices, then it may be that neither midpoint contacts them.
						//so in that case, we can still check if two non-equal vertices are contained in the line.
						if (p != (*this)[first_pass_vertex]) {
							//the line passes through two vertices.
							//so long as they do not form an edge of the hull, then this is an exclusive collision.
							
							Point<2, T> midpoint = ((*this)[first_pass_vertex] + p) / 2;
							if (this->contains(midpoint, polygon_out_this, inclusive)) {
								if (polygon_out_other) {
									polygon_out_other->clear();
									polygon_out_other->addVertex(line[0]);
									polygon_out_other->addVertex(line[1]);
								}
								return true;
							}
							
						}
					}
					//now form midpoints between this vertex and the two other ends of the line.
					if (line.a != p) {
						Point<2, T> midpoint = (line.a + p) / 2;
						if (this->contains(midpoint, polygon_out_this, inclusive)) {
							if (polygon_out_other) {
								polygon_out_other->clear();
								polygon_out_other->addVertex(line[0]);
								polygon_out_other->addVertex(line[1]);
							}
							return true;
						}
					}
					if (line.b != p) {
						Point<2, T> midpoint = (line.b + p) / 2;
						if (this->contains(midpoint, polygon_out_this, inclusive)) {
							if (polygon_out_other) {
								polygon_out_other->clear();
								polygon_out_other->addVertex(line[0]);
								polygon_out_other->addVertex(line[1]);
							}
							return true;
						}
					}
					
				}
			}
			
			//check to see if there are two vertices that are on the line segment (inclusively).
			//if so, then widen the containing polygon wtih vertices to the left and right of the line, if possible (if it isn't possible, then there is no collision).
			/*int first = -1;
			 int second =-1;
			 for(int i = 0;i<this->size();i++){
			 Point<2,T> & p = (*this)[i];
			 if (first<0 || (*this)[first]!=p){
			 if(line.contains(p,true)){//this check is intentionally inclusive
			 if(first<0){
			 first=i;
			 }else{
			 second = i;
			 break;
			 }
			 }
			 }
			 }
			 if(first>=0 && second>=0){
			 Point<2,T> & p1 =  (*this)[first];
			 Point<2,T> & p2 =  (*this)[second];

			 int outer = -1;
			 int inner = -1;
			 //then this MAY be a (non-inclusive) collision, so long as p1 p2 is not an edge of the polygon.
			 for(int i = first+1;i<second;i++){
			 Point<2,T> & p = (*this)[i];
			 if(line.whichSide(p)!=0){
			 inner = i;
			 break;
			 }
			 }
			 if(inner>-1){
			 for(int i = second+1; i< this->size()+first;i++){
			 Point<2,T> & p = (*this)[i];
			 if(line.whichSide(p)!=0){
			 outer = i;
			 break;
			 }
			 }

			 if(outer >-1){
			 //then the line passes through the convex polygon exactly at a vertex, but is not collinear with an edge.
			 //further more, these 4 points form an intersecting polygon.
			 r=true;
			 if(polygon_out_this){
			 polygon_out_this->clear();
			 polygon_out_this->addVertex((*this)[first]);
			 polygon_out_this->addVertex((*this)[inner]);
			 polygon_out_this->addVertex((*this)[second]);
			 polygon_out_this->addVertex((*this)[outer]);
			 }
			 if(polygon_out_other){
			 polygon_out_other->clear();
			 polygon_out_other->addVertex(line[0]);
			 polygon_out_other->addVertex(line[1]);
			 }

			 }

			 }

			 }*/
		}
		
#ifndef NDEBUG
		NConvexPolygon<2, T> test;
		test.addVertex(line.a);
		test.addVertex(line.b);
		if (this->intersects2d(test, nullptr, nullptr, inclusive) != r) {
			std::cout << line << "\n";
			std::cout << (*this) << "\n";
			edgesIntersectLine2d(line, nullptr, inclusive);
		}
		assert(this->intersects2d(test, nullptr, nullptr, inclusive) == r);
#endif
		return r;
		
	} else if (shape.getType() == CONVEX_POLYGON) {
		static int iterp = 0;
		if (++iterp == 9637275) {
			int a = 1;
		}
		
		ConvexPolygon<2, T> & c = (ConvexPolygon<2, T>&) shape;
		if (c.size() < this->size()) {
			return c.intersects(*this, polygon_out_other, polygon_out_this, inclusive);
		}
		
		if (this->size() == 0 || c.size() == 0) {
			return false;
		} else if (this->size() == 1) {
			if (c.contains((*this)[0], polygon_out_other, inclusive)) {
				if (polygon_out_this) {
					polygon_out_this->clear();
					polygon_out_this->addVertex((*this)[0]);
				}
				return true;
			} else {
				return false;
			}
			/*			if(!inclusive){
			 return false;
			 }
			 if(c.size()==1){
			 //then the two vertices are considered to collide if they are identical
			 return (*this)[0]==c[0];
			 }else{
			 return c.contains((*this)[0],inclusive);
			 }*/
		} else if (c.size() == 1) {
			if (!inclusive)
				return false;
			if (contains(c[0], polygon_out_this, inclusive)) {
				if (polygon_out_other) {
					polygon_out_other->clear();
					polygon_out_other->addVertex(c[0]);
				}
				return true;
			} else {
				return false;
			}
		}
		
		ConvexPolygon<2, T> & w = (ConvexPolygon<2, T>&) *this;
		
		//Separating Axis Theorem for collision detection between two convex polygons
		//loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
		//If any of the projections are non-intersection, then these don't collide; else, they do collide
		if (this->size() > 1) {
			for (int i = 0; i < this->size(); i++) {
				auto & cur = (*this)[i];
				auto & prev = (*this)[i - 1];
				Point<2, T> edge = cur - prev;
				Point<2, T> un_normalized_normal(-edge.y, edge.x);
				
				//now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
				//Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
				//(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
				T left = numeric<T>::infinity();
				T right = -numeric<T>::infinity();
				for (auto & p : *this) {
					T projection = un_normalized_normal.dot(p);
					if (projection < left) {
						left = projection;
					}
					if (projection > right) {
						right = projection;
					}
				}
				
				bool seenLeft = false;
				bool seenRight = false;
				for (auto & p : c) {
					T projection = un_normalized_normal.dot(p);
					if (inclusive) {
						if (projection >= left && projection <= right) {
							seenRight = true;
							seenLeft = true;
							break;
						} else if (projection < left) {
							seenLeft = true;
							if (seenRight) {
								break;
							}
						} else if (projection > right) {
							seenRight = true;
							if (seenLeft) {
								break;
							}
						} else if (seenLeft && projection > left) {
							seenRight = true;
							break;
						} else if (seenRight && projection < right) {
							seenRight = true;
							break;
						}
					} else {
						if (projection > left) {
							seenRight = true;
							if (seenLeft) {
								break;
							}
						}
						if (projection < right) {
							seenLeft = true;
							if (seenRight) {
								break;
							}
						}
						
					}
				}
				if (!(seenLeft && seenRight)) {
					return false;
				}
			}
		}
		if (c.size() > 1) {
			//now test the axis produced by the other polygon
			for (int j = 0; j < c.size(); j++) {
				auto & cur = c[j];
				auto & prev = c[j - 1];
				Point<2, T> edge = cur - prev;
				Point<2, T> un_normalized_normal(-edge.y, edge.x);
				
				T left = numeric<T>::infinity();
				T right = -numeric<T>::infinity();
				for (auto & p : c) {
					T projection = un_normalized_normal.dot(p);
					if (projection < left) {
						left = projection;
					}
					if (projection > right) {
						right = projection;
					}
				}
				bool seenLeft = false;
				bool seenRight = false;
				
				for (auto & p : *this) {
					T projection = un_normalized_normal.dot(p);
					if (inclusive) {
						if (projection >= left && projection <= right) {
							seenRight = true;
							seenLeft = true;
							break;
						} else if (projection < left) {
							seenLeft = true;
							if (seenRight) {
								break;
							}
						} else if (projection > right) {
							seenRight = true;
							if (seenLeft) {
								break;
							}
						} else if (seenLeft && projection > left) {
							seenRight = true;
							break;
						} else if (seenRight && projection < right) {
							seenRight = true;
							break;
						}
					} else {
						if (projection > left) {
							seenRight = true;
							if (seenLeft) {
								break;
							}
						}
						if (projection < right) {
							seenLeft = true;
							if (seenRight) {
								break;
							}
						}
						
					}
				}
				if (!(seenLeft && seenRight)) {
					return false;
				}
			}
			//If no axis overlapped, then they did in fact intersect
			
			//find minimal subsets of the polygons that are sufficient to intersect
			if (polygon_out_this || polygon_out_other) {
				if (polygon_out_this)
					polygon_out_this->clear();
				if (polygon_out_other)
					polygon_out_other->clear();
				
				//there are several possibilities for collisions, depending on whether we are inclusive of edges/vertices or not.
				ConvexPolygon<2, T> & h1 = *this;
				ConvexPolygon<2, T> & h2 = c;
				
				for (int i = 0; i < h1.size(); i++) {
					Point<2, T> & prev = h1[i - 1];
					Point<2, T> & p = h1[i];
					LineSegment<2, T> edge1(prev, p);
					for (int j = 0; j < h2.size(); j++) {
						Point<2, T> & prev2 = h2[j - 1];
						Point<2, T> & p2 = h2[j];
						LineSegment<2, T> edge2(prev2, p2);
						if ((inclusive || !edge1.collinear(edge2.a, edge2.b)) && edge1.intersects(edge2, inclusive)) {
							if (polygon_out_this) {
								polygon_out_this->addVertex(prev);
								polygon_out_this->addVertex(p);
								assert(polygon_out_this->orderClockwise());
							}
							if (polygon_out_other) {
								polygon_out_other->addVertex(prev2);
								polygon_out_other->addVertex(p2);
								assert(polygon_out_other->orderClockwise());
							}
							
							return true;
						}
					}
				}
				
				//if no intersecting line segment was found, then it follows that one of the polygons is wholly contained in the other (unless collisions are not inclusive, in which case an edge may be passing through one or more vertices).
				//so treat this as a contained point problem - pick one of the points in the contained polygon (arbitrarily), and a containing triangle
				//from the other polygon, and learn that one of these 4 points must be disabled.
				//Note that this may fail if collisions are not inclusive of the edges/vertices.
				for (int i = 0; i < h1.size(); i++) {
					if (h2.contains(h1[i], polygon_out_other, inclusive)) {
						if (polygon_out_this) {
							polygon_out_this->addVertex(h1[i]);
							assert(polygon_out_this->orderClockwise());
						}
						return true;
					}
				}
				
				for (int i = 0; i < h2.size(); i++) {
					if (h1.contains(h2[i], polygon_out_this, inclusive)) {
						if (polygon_out_other) {
							polygon_out_other->addVertex(h2[i]);
							assert(polygon_out_other->orderClockwise());
						}
						return true;
					}
				}
				
				assert(!inclusive);
				//if collisions are not inclusive, then there are three more possibilities:
				//1) one of the polygons has an edge that cleaves the other polygon in two, but whose endpoints are exactly on the edges of the other polygon
				//2) one of the polygons is a line (or has only two distinct points), and that line passes through or ends on two vertices of the other polygon.
				//3) the two polygons are identical. that possibility can be accounted for by finding any three identical points in the two polygons.
				
				//1) find that cleaving edge, and then find, from the other polygon, two vertices that form a line that crosses through the first line (can we really always find two such vertices?).
				//
				for (int hull = 0; hull < 2; hull++) {
					ConvexPolygon<2, T> & hull1 = hull ? h2 : h1;
					ConvexPolygon<2, T> & hull2 = hull ? h1 : h2;
					NConvexPolygon<2, T> * polygon_out_1 = hull ? polygon_out_other : polygon_out_this;
					NConvexPolygon<2, T> * polygon_out_2 = hull ? polygon_out_this : polygon_out_other;
					LineSegment<2, T> testline;
					for (int i = 0; i < hull1.size(); i++) {
						
						testline.a = hull1[i - 1];
						testline.b = hull1[i];
						
						bool intersect = hull2.intersects(testline, polygon_out_2, polygon_out_1, inclusive);
						if (intersect) {
							return true;
						}
						/*	Point<2,T> mid = (hull1[i-1] +  hull1[i])/2;//testing the midpoint is only guaranteed to work if the line has endpoints exactly on the other hull.


						 if(hull2.contains(mid,inclusive)){
						 //then this is may be a cleaving line.
						 if(polygon_out_1){
						 polygon_out_1->clear();
						 polygon_out_1->addVertex(hull1[i-1]);
						 polygon_out_1->addVertex(hull1[i]);

						 }
						 //now find two vertices of h2 that form a line that intersects this one. the vertices will come from separate sides of h2.
						 if(polygon_out_2){
						 polygon_out_2->clear();
						 static LineSegment<2,T> line1;

						 line1.a = hull1[i-1];
						 line1.b = hull1[i];
						 bool found_left=false;
						 bool found_right=false;
						 for (int j = 0;j<hull2.size();j++){
						 Point<2,T> & p = hull2[j];
						 int side = line1.whichSide(p);
						 if(side<0 && ! found_left){
						 found_left=true;
						 polygon_out_2->addVertex(p);
						 if(found_right){
						 break;
						 }
						 }else if(side>0 && ! found_right){
						 found_right=true;
						 polygon_out_2->addVertex(p);
						 if(found_left){
						 break;
						 }
						 }
						 }
						 assert(polygon_out_2->orderClockwise());
						 if(found_left && found_right){
						 assert(line1.intersects(*polygon_out_2,true));
						 //there is still one more possibility - which is that the line segment we have found through the vertices to the left and right of the cleaving line
						 //cross the cleaving line exactly at an endpoint.
						 //in that case, if the collision is not inclusive, then we can resolve this by finding any point from h2 that is not collinear with that line.
						 if(!inclusive && !line1.intersects(*polygon_out_2,false)){
						 static LineSegment<2,T> line2;
						 line2.a = (*polygon_out_2)[0];
						 line2.b = (*polygon_out_2)[1];
						 assert(hull2.dbg_orderClockwise());
						 assert(hull2.orderClockwise());
						 bool found=false;
						 for (int j = 0;j<hull2.size();j++){
						 Point<2,T> & p = hull2[j];
						 int side = line2.whichSide(p);
						 if(side!=0){
						 polygon_out_2->addVertex(p);
						 //We need to reorder the vertices of polygon 2 to put them back into clockwise order if we add a point.
						 polygon_out_2->reorderVertices();
						 found=true;
						 break;
						 }
						 }

						 assert(polygon_out_2->orderClockwise());
						 //assert(polygon_out_2->dbg_orderClockwise());
						 //polygon_out_2->reorderVertices();
						 assert(polygon_out_2->dbg_orderClockwise());
						 assert(found);
						 assert(polygon_out_2->intersects(line1,false));
						 if(!found){
						 exit(4);//this shouldn't be possible for convex polygons!
						 }
						 }
						 }else{
						 //then we need to look for a different line.
						 assert(false);//is this actually a reachable case? It shouldn't be, because we are only reaching this code if polygon_out_other is non-null, and that is an optional parameter.
						 exit(4);
						 }
						 return true;
						 }
						 }*/
					}
				}
				
				if (polygon_out_this) {
					polygon_out_this->clear();
				}
				if (polygon_out_other) {
					polygon_out_other->clear();
				}
				
				//one (hopefully) final possibility is that
				//3) the two polygons are identical. pick any three identical points to form two intersecting triangles.
				int count = 0;
				//this probably doesn't have to be quadratic
				for (int i = 0; i < h1.size(); i++) {
					Point<2, T> & p = h1[i];
					
					for (int j = 0; j < h2.size(); j++) {
						Point<2, T> & p2 = h2[j];
						if (p == p2) {
							if (polygon_out_this) {
								polygon_out_this->addVertex(h1[i]);
							}
							if (polygon_out_other) {
								polygon_out_other->addVertex(h2[j]);
							}
							if (++count == 3) {
								polygon_out_other->reorderVertices();
								polygon_out_this->reorderVertices();
								return true;
							}
							break;
						}
					}
				}
				if (polygon_out_this) {
					polygon_out_this->clear();
				}
				if (polygon_out_other) {
					polygon_out_other->clear();
				}
				
				std::cout << "Failed to find intersecting polygon for intersection " << iterp << ", aborting.\n";
				if (inclusive)
					std::cout << "(inclusive)\n";
				if (ignore_vertices)
					std::cout << "(ignore_vertices)\n";
				std::cout << h1 << "\n";
				std::cout << h2 << "\n";
				
				assert(false);
				exit(4);
			}
			
			return true;
		}
	}
	assert(false);
	exit(4);
	return false;
}

//put the vertices into clockwise order
template<unsigned int D, class T>
void NConvexPolygon<D, T>::reorderVertices2d() {
	this->vertices_clockwise = true;
	if (vertices.size() <= 2) {
		return;
	}
	T centerX = 0;
	T centerY = 0;
	std::vector<Point<2, T>> & w = (std::vector<Point<2, T>> &) vertices;
	for (auto & p : w) {
		centerX += p.x;
		centerY += p.y;
	}
	T npoints =  (int)vertices.size();
	centerX /=npoints;
	centerY /=npoints;
	
	//from http://stackoverflow.com/a/6989383
	
	struct clockwise_lt {
		const std::vector<Point<2, T>> & points;
		T centerX;
		T centerY;
		clockwise_lt(const std::vector<Point<2, T>> & points, T centerX, T centerY) :
				points(points), centerX(centerX), centerY(centerY) {
			
		}
		
		bool operator ()(int id_a, int id_b) {
			assert(id_a >= 0);
			assert(id_a < points.size());
			assert(id_b >= 0);
			assert(id_b < points.size());
			auto & a = points[id_a];
			auto & b = points[id_b];
			if (a[0] - centerX >= 0 && b[0] - centerX < 0)
				return true;
			if (a.x - centerX < 0 && b.x - centerX >= 0)
				return false;
			if (a.x - centerX == 0 && b.x - centerX == 0) {
				if (a.y - centerY >= 0 || b.y - centerY >= 0)
					return a.y > b.y;
				return b.y > a.y;
			}
			
			// compute the cross product of vectors (center -> a) x (center -> b)
			T det = (a.x - centerX) * (b.y - centerY) - (b.x - centerX) * (a.y - centerY);
			if (det < 0)
				return true;
			if (det > 0)
				return false;
			
			// points a and b are on the same line from the center
			// check which point is closer to the center
			T d1 = (a.x - centerX) * (a.x - centerX) + (a.y - centerY) * (a.y - centerY);
			T d2 = (b.x - centerX) * (b.x - centerX) + (b.y - centerY) * (b.y - centerY);
			return d1 > d2;
		}
	};
	//this should ideally be avoided...
	static std::vector<int> points_clockwise;
	points_clockwise.clear();
	for (int i = 0; i < vertices.size(); i++) {
		points_clockwise.push_back(i);
	}
	std::sort(points_clockwise.begin(), points_clockwise.end(), clockwise_lt(vertices, centerX, centerY));
	//do this in place later
	static std::vector<Point<2, T>> oldPoints;
	oldPoints.clear();
	for (int i = 0; i < vertices.size(); i++) {
		oldPoints.push_back(vertices[i]);
	}
	for (int i = 0; i < vertices.size(); i++) {
		vertices[i] = oldPoints[points_clockwise[i]];
	}

	assert(this->dbg_orderClockwise());
}

#endif /* CONVEXPOLYGON_H_ */
