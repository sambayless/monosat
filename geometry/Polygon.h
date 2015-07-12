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

#ifndef POLYGON_H_
#define POLYGON_H_
#include <math.h>
#include "Shape.h"
#include <gmpxx.h>
#include "bounds/BoundingVolume.h"
#include <iostream>
/**
 * A concrete polygon (or, for D>2, a polytope)
 */
template<unsigned int D, class T>
class Polygon: public Shape<D, T> {
public:
	//List of vertices in clockwise order
	
	BoundingVolume<D, T> * bound = nullptr;
protected:
	int distinct_vertices = -1;
	bool vertices_clockwise = false;
	bool bounds_uptodate = false;
public:
	Polygon() {
		
	}
	
	virtual ~Polygon() {
		if (bound != nullptr) {
			delete (bound);
			bound = nullptr;
		}
	}
	;
	virtual ShapeType getType() {
		return POLYGON;
	}
	void setBoundingVolume(BoundingVolume<D, T> * bound) {
		this->bound = bound;
	}
	
	bool hasBound() const {
		return bound;
	}
	
	BoundingVolume<D, T> * getBound() {
		if (!bounds_uptodate) {
			bounds_uptodate = true;
			bound->update();
		}
		return bound;
	}
	
	virtual bool contains(const Point<D, T> & point, bool inclusive = true) {
		if (!boundContains(point, inclusive)) {
			return false;
		}
		if (D == 2) {
			return contains2d((const Point<2, T> &) point, inclusive);
		}
		assert(false);
		return false;
	}
	
	//virtual bool findContainingConvex(const Point<D,T> & point,NPolygon<D,T> & polygon_out)=0;
	
	virtual bool intersects(Shape<D, T> & s, bool inclusive = true) {
		assert(false);
		return false;
	}
	
	virtual int size() const=0;
	virtual void update()=0;

	//Returns the vertices of the polygon, in clockwise order.
	virtual const Point<D, T>& operator [](int index) const =0;
	virtual Point<D, T>& operator [](int index)=0;

	virtual Point<D, T>& back() {
		return (*this)[size() - 1];
	}
	
	virtual T getArea() {
		if (D == 2) {
			return getArea2d();
		} else
			assert(false);
		return 0;
	}
	virtual T getPerimeter() {
		if (D == 2) {
			return getPerimeter2d();
		} else
			assert(false);
		return 0;
	}
	
	int distinctVertices() {
		if (distinct_vertices < 0) {
			if (!vertices_clockwise) {
				int same_count = 0;
				for (int i = 0; i < size(); i++) {
					for (int j = i + 1; j < size(); j++) {
						if ((*this)[i] == (*this)[j]) {
							same_count++;
						}
					}
				}
				distinct_vertices = size() - same_count;
			} else {
				assert(vertices_clockwise);
				assert(dbg_orderClockwise());
				distinct_vertices = 0;
				for (int i = 0; i < size(); i++) {
					if ((*this)[i] != (*this)[i - 1]) {
						distinct_vertices++;
					}
				}
			}
		}
		assert(dbg_distinctVertexCount());
		return distinct_vertices;
	}
	
	inline bool boundContains(const Point<D, T> & p, bool inclusive = true) {
		if (!bound)
			return true;
		if (!bounds_uptodate) {
			bounds_uptodate = true;
			bound->update();
		}
		return bound->contains(p, inclusive);
	}
	inline bool boundIntersects(Shape<D, T> & p, bool inclusive = true) {
		if (!bound)
			return true;
		if (!bounds_uptodate) {
			bounds_uptodate = true;
			bound->update();
		}
		return bound->intersects(p, inclusive);
	}
	static bool pointInTriangle2d(const Point<D, T> & p, const Point<D, T> & p0, const Point<D, T> & p1,
			const Point<D, T> &p2, bool inclusive) {
		//from http://stackoverflow.com/a/14382692
		assert(dbg_orderClockwise2dTri(p0, p1, p2));
		
		T s = (p2.y * p0.x - p2.x * p0.y + (p0.y - p2.y) * p.x + (p2.x - p0.x) * p.y);
		T t = (p2.x * p1.y - p2.y * p1.x + (p2.y - p1.y) * p.x + (p1.x - p2.x) * p.y);
		T area2 = (-p1.y * p0.x + p2.y * (-p1.x + p0.x) + p2.x * (p1.y - p0.y) + p1.x * p0.y);
		assert(area2 >= 0);
		if (inclusive) {
			return s >= 0 && t >= 0 && (s + t <= area2);
		} else {
			return s > 0 && t > 0 && (s + t < area2);
		}
		//from http://stackoverflow.com/a/13301035
		/*	T dif23y = p2.y - p3.y;
		 T dif32x = p3.x - p2.x;
		 T dif13y = p1.y - p3.y;
		 T dif13x = p1.x - p3.x;
		 T dif03x = p.x - p3.x;
		 T dif03y = p.y - p3.y;
		 T alpha = (dif23y*(dif03x) + dif32x*(dif03y)) /
		 (dif23y*(dif13x) + dif32x*(dif13y));
		 T beta = ((p3.y - p1.y)*(dif03x) + (dif13x)*(dif03y)) /
		 (dif23y*(dif13x) + dif32x*(dif13y));
		 //T gamma = 1 - alpha - beta;
		 return (alpha>=0 && beta>=0 && (( alpha + beta)<=1) ) //&& gamma>=0);//should these be > or >=? for inclusive or exclusive containment?*/

	}
	
	bool orderClockwise() const {
		if (D == 2) {
			return orderClockwise2d();
		}
		assert(false);
		return true;
	}
	bool orderClockwise2d() const;
protected:
	bool dbg_orderClockwise() {
		if (D == 2) {
			return dbg_orderClockwise2d();
		}
		return true;
	}
	bool dbg_orderClockwise2d() const;
	bool dbg_boundsUpToDate() {
		return true;
	}
	
	bool dbg_distinctVertexCount() const {
#ifndef NDEBUG
		if (distinct_vertices >= 0) {
			int same_count = 0;
			for (int i = 0; i < size(); i++) {
				for (int j = i + 1; j < size(); j++) {
					if ((*this)[i] == (*this)[j]) {
						same_count++;
					}
				}
			}
			assert(size() - same_count == distinct_vertices);
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
protected:
	
	T getArea2d();
	T getPerimeter2d();
	//put the vertices into clockwise order
	void reorderVertices2d();

	bool contains2d(const Point<2, T> & point, bool inclusive);

	//Note: for convenience, the first point of the wrap is also the last point (it is duplicated).
	
	struct wrap_iterator {
		Point<2, T> &operator*() {
			return verts[pos % verts.size()];
		}
		
		wrap_iterator operator++() {
			wrap_iterator i = *this;
			pos++;
			return i;
		}
		wrap_iterator operator++(int ignore) {
			pos++;
			return *this;
		}
		Point<2, T>* operator->() {
			return &verts[pos % verts.size()];
		}
		bool operator==(const wrap_iterator& rhs) {
			return pos == rhs.pos;
		}
		bool operator!=(const wrap_iterator& rhs) {
			return pos != rhs.pos;
		}
		
		Polygon<D, T> & verts;
		int pos = 0;

		wrap_iterator(Polygon<D, T> & verts, int pos) :
				verts(verts), pos(pos) {
		}
		;
		// private constructor for begin, end
	};

	typedef wrap_iterator iterator;

public:
	iterator begin() {
		return iterator(*this, 0);
	}
	
	iterator end() {
		return iterator(*this, size());
	}
	
	iterator end_wrap() {
		return iterator(*this, size() + 1);
	}
	
};
template<unsigned int D, class T>
class NPolygon: public Polygon<D, T> {
public:
	//List of vertices in clockwise order
	
	std::vector<Point<D, T>> vertices;

	NPolygon() {
		
	}
	
	virtual ~NPolygon() {
		
	}
	explicit NPolygon(const NPolygon&from) {
		vertices = from.vertices;
	}
	
	virtual ShapeType getType() {
		return POLYGON;
	}
	
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
			if (p != vertices[0] && p != vertices.last()) {
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
		vertices.pop();
		this->bounds_uptodate = false;
	}
	
	void clearVertices() {
		vertices.clear();
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
	
	virtual bool isConvex() {
		bool seenPositive = false;
		bool seenNegative = false;
		for (int i = 0; i < size(); i++) {
			const Point<2, T> & prev = i > 0 ? (*this)[i - 1] : (*this).back();
			const Point<2, T> & p = (*this)[i];
			const Point<2, T> & next = i < size() - 1 ? (*this)[i + 1] : (*this)[0];
			Point<2, T> a = p - prev;
			Point<2, T> b = next - p;
			T s = cross2d(a, b);
			seenPositive |= s > 0;
			seenNegative |= s < 0;
			if (seenPositive && seenNegative)
				return false;
		}
		return true;
	}
private:
	// copy ops are private to prevent copying
	//Polygon(const Polygon& from); // no implementation
	NPolygon& operator=(const NPolygon& from); // no implementation
			
	//put the vertices into clockwise order
	void reorderVertices2d();
	
};
template<unsigned int D, class T>
std::ostream & operator<<(std::ostream & str, NPolygon<D, T> & polygon) {
	str << "Polygon=[";
	for (const auto & p : polygon) {
		str << p << ",";
	}
	str << "]";
	return str;
}
template<unsigned int D, class T>
std::ostream & operator<<(std::ostream & str, Polygon<D, T> & polygon) {
	str << "Polygon=[";
	for (const auto & p : polygon) {
		str << p << ",";
	}
	str << "]";
	return str;
}
template<unsigned int D, class T>
inline bool Polygon<D, T>::orderClockwise2d() const {
	
	//from http://stackoverflow.com/a/1165943
	//but reversed relative to stackoverflow, to work for clockwise instead of counterclockwise
	T sum = 0;
	for (int i = 0; i < size(); i++) {
		Point<2, T> & a = (Point<2, T> &) (*this)[i - 1]; //(i>0? (*this)[i-1]:(*this)[this->size()-1]);
		Point<2, T> & b = (Point<2, T> &) (*this)[i];
		sum += (b.x - a.x) * (b.y + a.y);
	}
	return sum >= 0;
}

template<unsigned int D, class T>
inline bool Polygon<D, T>::dbg_orderClockwise2d() const {
#ifndef NDEBUG
	//from http://stackoverflow.com/a/1165943
	if (vertices_clockwise) {
		T sum = 0;
		for (int i = 0; i < size(); i++) {
			Point<2, T> & a = (Point<2, T> &) (*this)[i - 1]; //(i>0? (*this)[i-1]:(*this)[this->size()-1]);
			Point<2, T> & b = (Point<2, T> &) (*this)[i];
			sum += (b.x - a.x) * (b.y + a.y);
		}
		assert(sum >= 0);
		if (this->size() == 3) {
			assert(dbg_orderClockwise2dTri((*this)[0], (*this)[1], (*this)[2]));
		}
	}
	
#endif
	return true;
}

template<unsigned int D, class T>
T Polygon<D, T>::getArea2d() {
	//std::vector<Point<2,T>> &  points = getVertices();
	int sz = size();
	T sum = 0;
	for (int i = 0; i < sz; i++) {
		Point<2, T>& prev = i > 0 ? (*this)[i - 1] : (*this)[sz - 1];
		Point<2, T>& cur = (*this)[i];
		sum += prev[0] * cur[1] - cur[0] * prev[1];
	}
	using namespace std;
	return abs(sum / 2.0);
}

template<unsigned int D, class T>
bool Polygon<D, T>::contains2d(const Point<2, T> & point, bool inclusive) {
//This is the PNPOLY test, adapted from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
//This method has the following copyright.
	/*	Copyright (c) 1970-2003, Wm. Randolph Franklin

	 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

	 Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
	 Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
	 The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission.

	 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
	 */
	int i;
	int j;
	bool result = false;
	for (i = 0, j = (*this).size() - 1; i < (*this).size(); j = i++) {
		if (((*this)[i].y > point.y) != ((*this)[j].y > point.y)
				&& (point.x
						< ((*this)[j].x - (*this)[i].x) * (point.y - (*this)[i].y) / ((*this)[j].y - (*this)[i].y)
								+ (*this)[i].x)) {
			result = !result;
		}
	}
	return result;
}
/**
 *from http://objectmix.com/graphics/314163-polygon-approximation-simplification-inner-outer-bounds.html
 two tricks:
 Segment trees -- Assuming that your point-in-polygon test is
 counting polygon segments that intersect a vertical ray eminating
 from that point, for each polygon, you sort the vertices by x-value
 and construct a tree in which each node represents a range of
 adjacent x values. Off of each tree node you then build a linked
 list of pointers to the polygon segments that straddle that
 interval. Thus, your point-in-polygon test need only traverse the
 tree to find the interval containing that point's x-value, and then
 loop through the short list of segments that lie vertically above
 or below that point. My testing shows (as one would expect)
 logrithmic performance improvement: a doubling of processing speed
 for polygons with 160 vertices and a 10x improvement at 1000
 vertices. I also determined that for polygons with fewer than 20
 vertices, the performance improvement didn't warrant the time spent
 constructing the segment tree, but of course, that's strictly a
 function of the data I'm processing. The set-up time becomes
 comparativley less significant as you process more points within
 that polygon.

 Gridding off the space -- Similar to the R-tree you've already
 implemented, but simpler and faster. You simply define a grid that
 partitions the range of x-y values, and then for each grid square,
 you build a linked list of pointers to the polygons that contact
 that square. If you make the grid small enough such that typical
 polygons will contain several grid squares, then you can even
 eliminate the need to run your point-in-polygon test (since every
 point within that "enclosed" grid will be within that polygon).
 This dramatically improves performance (by my testing of ad valorem
 tax data, nearly two orders of magnitude) over the brute-force
 method of testing every point against every polygon. As a further
 enhancement, if you make your grid size a power of two and if your
 test points and polygon vertices are specified by integer values,
 then you might even use the C-language shift operator (thus
 avoiding division) to find the x-y index values of the grid
 containing that point. The price you pay for this point-in-polygon
 performance gain is memory and set up time. If the polygon data is
 relatively static (as are the city limits, etc. I process), then
 you can archive the grid-to-polygon lists and directly load that
 into memory on subsequence runs of new points against those same
 polygons.

 The data you're processing will greatly influence performance. For
 a relatively small number of complex polygons processed against
 millions of points, the two tricks I've proposed yield very
 dramatic performance improvement. If your processing involves a
 large number of small polyons (polygons with a small number of
 vertices) being run on a small number of points (you mentioned
 "thousands"), then you probably shouldn't bother constructing
 segment trees (except perhaps for the largest of polygons).
 */

//Note that this is subject to rounding errors.
template<unsigned int D, class T>
T Polygon<D, T>::getPerimeter2d() {
	
	T sum = 0;
	for (int i = 1; i < size(); i++) {
		Point<2, T> prev = (*this)[i - 1];
		Point<2, T> cur = (*this)[i];
		T xdist = cur[0] - prev[0];
		T ydist = cur[1] - prev[1];
		T v = xdist * xdist + ydist * ydist;
		T s = sqrt(v);
		sum += s;
	}
	return sum;
}

//put the vertices into clockwise order
template<unsigned int D, class T>
void NPolygon<D, T>::reorderVertices2d() {
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
	centerX /= vertices.size();
	centerY /= vertices.size();
	
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
	assert(vertices.size() > 0);
	
	for (int i = 0; i < vertices.size(); i++) {
		oldPoints.push_back(vertices[i]);
	}
	for (int i = 0; i < vertices.size(); i++) {
		vertices[i] = oldPoints[points_clockwise[i]];
	}
	
	assert(this->dbg_orderClockwise());
	assert(this->dbg_distinctVertexCount());
}

#endif /* POLYGON_H_ */
