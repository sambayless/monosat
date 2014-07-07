/*
 * Polygon.h
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#ifndef POLYGON_H_
#define POLYGON_H_
#include <math.h>
#include "Shape.h"
#include <gmpxx.h>
#include "bounds/BoundingVolume.h"
/**
 * A concrete polygon (or, for D>2, a polytope)
 */
template<unsigned int D,class T>
class Polygon:public Shape<D,T>{
public:
	//List of vertices in clockwise order

	std::vector<Point<D,T>> vertices;

	BoundingVolume<D,T> * bound=nullptr;
	bool vertices_clockwise=false;
	bool bounds_uptodate=false;
	Polygon(){

	}
	explicit Polygon(const Polygon&from){
		vertices=from.vertices;
	}

	virtual ~Polygon(){
		if(bound!=nullptr){
			delete(bound);
			bound=nullptr;
		}
	};
	virtual ShapeType getType(){
		return POLYGON;
	}
	void setBoundingVolume(BoundingVolume<D,T> * bound){
		this->bound=bound;
	}


	virtual bool contains(const Point<D,T> & point){
		if(!boundContains(point)){
			return false;
		}
		if(D==2){
			return contains2d((const Point<2,T> &)point);
		}
		assert(false);
		return false;
	}

	virtual bool findContainingConvex(const Point<D,T> & point,Polygon<D,T> & polygon_out){
		assert(false);
		return false;
	}

	virtual bool intersects(Shape<D,T> & s){
		assert(false);
		return false;
	}

	int size()const {
		return vertices.size();
	}

	void update(){
		if(!vertices_clockwise){
			reorderVertices();
		}
		bounds_uptodate=false;
	}

	void clear(){
		vertices.clear();
	}

	int size(){
		return vertices.size();
	}
	void addVertex(Point<D,T> p){
		vertices_clockwise=false;
		bounds_uptodate=false;
		vertices.push_back(p);
	}
/*	void addVertex(Point<D,T> & p){
		vertices_clockwise=false;
		bounds_uptodate=false;
		vertices.push_back(p);
	}*/
	//add a vertex, assuming that it will preserve clockwise order
	void addVertexUnchecked(Point<D,T>  p){
		vertices.push_back(p);
		assert(dbg_orderClockwise());
		assert(dbg_boundsUpToDate());
	}

	void popVertex(){
		vertices.pop();
		bounds_uptodate=false;
	}

	void clearVertices(){
		vertices.clear();
		bounds_uptodate=false;
	}

	//Returns the vertices of the polygon, in clockwise order.

	std::vector<Point<D,T> > & getVertices(){
		if(!vertices_clockwise){
			reorderVertices();
		}
		dbg_orderClockwise();
		return vertices;
	}
    const Point<D,T>& operator [] (int index) const {
    	index = index %vertices.size();
    	if(index<0){
    		index+=vertices.size();
    	}
    	assert(index>=0);assert(index<vertices.size());
    	return vertices[index];
    }
    Point<D,T>&       operator [] (int index)       {
    	index = index %vertices.size();
    	if(index<0){
    		index+=vertices.size();
    	}
    	assert(index>=0);assert(index<vertices.size());
    	return vertices[index];
    }

	virtual T getArea(){
		if(D==2){
			return getArea2d();
		}else
			assert(false);
		return 0;
	}
	virtual T getPerimeter(){
		if(D==2){
			return getPerimeter2d();
		}else
			assert(false);
		return 0;
	}
	//put the vertices into clockwise order
	void reorderVertices(){
		if(D==2){
			reorderVertices2d();
		}else
			assert(false);
	}

	inline bool boundContains(const Point<D,T> & p){
		if(!bound)
			return true;
		if(!bounds_uptodate){
			bounds_uptodate=true;
			bound->update();
		}
		return bound->contains(p);
	}
	static bool pointInTriangle2d(const Point<D,T> & p,const Point<D,T> & p0,const Point<D,T> & p1, const Point<D,T> &p2){
		//from http://stackoverflow.com/a/14382692
		assert(dbg_orderClockwise2dTri(p0,p1,p2));

		T s = (p2.y*p0.x - p2.x*p0.y + (p0.y - p2.y)*p.x + (p2.x - p0.x)*p.y);
		T t = (p2.x*p1.y - p2.y*p1.x + (p2.y - p1.y)*p.x + (p1.x - p2.x)*p.y);
		T area2 =  (-p1.y*p0.x + p2.y*(-p1.x + p0.x) + p2.x*(p1.y - p0.y) + p1.x*p0.y);
		assert(area2>=0);
		return s>=0 && t>=0 && (s+t<=area2);
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
protected:
	bool dbg_orderClockwise(){
		return true;
	}

	bool dbg_boundsUpToDate(){
		return true;
	}

	static bool dbg_orderClockwise2dTri(Point<2,T> p1,Point<2,T> p2,Point<2,T> p3){
	#ifndef NDEBUG
		std::vector<Point<2,T>> points;
		points.push_back(p1);
		points.push_back(p2);
		points.push_back(p3);


			T sum = 0;
			for(int i = 0;i<points.size();i++){
				Point<2,T> & a = i>0? points[i-1]:points.back();
				Point<2,T> & b = points[i];
				sum+= (b.x - a.x)*(b.y+a.y);
			}
			assert(sum>=0);
			return sum>=0;

	#endif
		return true;
	}
private:
	// copy ops are private to prevent copying
	//Polygon(const Polygon& from); // no implementation
	Polygon& operator=(const Polygon& from); // no implementation
	T getArea2d();
	T getPerimeter2d();
	//put the vertices into clockwise order
	void reorderVertices2d();

	bool contains2d(const Point<2,T> & point);

	//Note: for convenience, the first point of the wrap is also the last point (it is duplicated).
	template< class ValueType >
	struct wrap_iterator {
	    ValueType &operator*() { return verts[pos%verts.size()]; }



	    wrap_iterator operator++() { wrap_iterator i = *this; pos++; return i; }
	    wrap_iterator operator++(int ignore) { pos++; return *this; }
	    ValueType* operator->() { return &verts[pos%verts.size()]; }
	    bool operator==(const wrap_iterator& rhs) { return pos == rhs.pos; }
	    bool operator!=(const wrap_iterator& rhs) { return pos != rhs.pos; }

	    int pos=0;
	    std::vector<ValueType> & verts;


	    wrap_iterator( std::vector<ValueType> & verts, int pos):verts(verts),pos(pos){}; // private constructor for begin, end
	};

	typedef wrap_iterator<  Point<D,T> > iterator;
	typedef wrap_iterator<  Point<D,T> const > const_iterator;
public:
	 iterator begin()
	{
		 return iterator(getVertices(),0);
	}

	iterator end()
	{
		return iterator(getVertices(),getVertices().size());
	}

	iterator end_wrap()
	{
		return iterator(getVertices(),getVertices().size()+1);
	}

/*	 const_iterator begin() const
	{
		 return const_iterator(vertices,0);
	}

	const_iterator end() const
	{
		return const_iterator(vertices,vertices.size());
	}

	const_iterator end_wrap() const
	{
		return const_iterator(vertices,vertices.size()+1);
	}*/
};

template<>
inline bool Polygon<2,double>::dbg_orderClockwise(){
#ifndef NDEBUG
	//from http://stackoverflow.com/a/1165943
	if(vertices_clockwise){
		double sum = 0;
		for(int i = 0;i<vertices.size();i++){
			Point<2,double> & a = i>0? vertices[i-1]:vertices.back();
			Point<2,double> & b = vertices[i];
			sum+= (b.x - a.x)*(b.y+a.y);
		}
		assert(sum>=0);
	}

#endif
	return true;
}

template<unsigned int D,class T>
T Polygon<D,T>::getArea2d(){
	std::vector<Point<2,T>> &  points = getVertices();

	T sum = 0;
	for (int i = 0;i<points.size();i++){
		Point<2,T>& prev = i>0? points[i-1]: points.back();
		Point<2,T>& cur = points[i];
		sum += prev[0]*cur[1]-cur[0]*prev[1];
	}
	return abs(sum/2.0);
}

template<unsigned int D,class T>
bool Polygon<D,T>::contains2d(const Point<2,T> & point){

std::vector< Point<2,T>> &  points =(std::vector< Point<2,T>> &) getVertices();
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
template<unsigned int D,class T>
T Polygon<D,T>::getPerimeter2d(){
	std::vector< Point<2,T>> &  w =(std::vector< Point<2,T>> &) getVertices();
	T sum = 0;
	for (int i = 1;i<w.size();i++){
		Point<2,T> prev = w[i-1];
		Point<2,T> cur = w[i];
		T xdist = cur[0]-prev[0];
		T ydist=cur[1]-prev[1];
		sum += sqrt(xdist*xdist + ydist*ydist);
	}
	return sum;
}

template<>
inline  double Polygon<2,double>::getPerimeter2d(){
	std::vector< Point<2,double>> &  w =(std::vector< Point<2,double>> &) getVertices();
	double sum = 0;
	for (int i = 1;i<w.size();i++){
		auto & prev = w[i-1];
		auto &  cur = w[i];
		double xdist = cur[0]-prev[0];
		double ydist=cur[1]-prev[1];
		sum += sqrt(xdist*xdist + ydist*ydist);
	}
	return sum;
}
template<>
inline mpq_class Polygon<2,mpq_class>::getPerimeter2d(){
	std::vector< Point<2,mpq_class>> &  w =(std::vector< Point<2,mpq_class>> &) getVertices();
	mpq_class sum = 0;
	for (int i = 1;i<w.size();i++){
		auto &  prev = w[i-1];
		auto &  cur = w[i];
		mpq_class xdist = cur[0]-prev[0];
		mpq_class ydist=cur[1]-prev[1];
		sum += sqrt((mpf_class)( xdist*xdist + ydist*ydist));
	}
	return sum;
}
//put the vertices into clockwise order
template<unsigned int D,class T>
void Polygon<D,T>::reorderVertices2d(){
	vertices_clockwise=true;
	if (vertices.size()<=2){
		return;
	}
	T centerX = 0;
	T centerY = 0;
	std::vector< Point<2,T>> &  w =(std::vector< Point<2,T>> &) vertices;
	for(auto & p:w){
		centerX+=p.x;
		centerY+=p.y;
	}
	centerX/=vertices.size();
	centerY/=vertices.size();

	//from http://stackoverflow.com/a/6989383

	struct clockwise_lt{
		const std::vector<Point<2,T>> & points;
		T centerX;
		T centerY;
		clockwise_lt(const std::vector<Point<2,T>> & points,T centerX, T centerY):points(points),centerX(centerX),centerY(centerY){

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
			T det = (a.x -centerX) * (b.y - centerY) - (b.x -centerX) * (a.y - centerY);
			if (det < 0)
				return true;
			if (det > 0)
				return false;

			// points a and b are on the same line from the center
			// check which point is closer to the center
			T d1 = (a.x -centerX) * (a.x -centerX) + (a.y - centerY) * (a.y - centerY);
			T d2 = (b.x -centerX) * (b.x -centerX) + (b.y - centerY) * (b.y - centerY);
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
	static std::vector<Point<2,T>> oldPoints;
	oldPoints.clear();
	for(int i =0;i<vertices.size();i++){
		oldPoints.push_back(vertices[i]);
	}
	for(int i =0;i<vertices.size();i++){
		vertices[i] = oldPoints[points_clockwise[i]];
	}


	assert(dbg_orderClockwise());
}


template<>
inline bool Polygon<2,mpq_class>::dbg_orderClockwise(){
#ifndef NDEBUG
	//from http://stackoverflow.com/a/1165943
	if(vertices_clockwise){
		mpq_class sum = 0;
		for(int i = 0;i<vertices.size();i++){
			Point<2,mpq_class> & a = i>0? vertices[i-1]:vertices.back();
			Point<2,mpq_class> & b = vertices[i];
			sum+= (b.x - a.x)*(b.y+a.y);
		}
		assert(sum>=0);
	}

#endif
	return true;
}





#endif /* POLYGON_H_ */
