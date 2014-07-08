/*
 * ConvexPolygon.h
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#ifndef CONVEXPOLYGON_H_
#define CONVEXPOLYGON_H_
#include "Polygon.h"
#include <gmpxx.h>
#include "core/Config.h"
#include "Line.h"
#include "LineSegment.h"

template<unsigned int D,class T>
class NConvexPolygon;

template<unsigned int D,class T>
class Triangle;

/**
 * A concrete, convex polygon (or, for D>2, a polytope)
 */
template<unsigned int D,class T>
class ConvexPolygon:public Polygon<D,T>{
public:
	static long stats_triangle_avoided;
	static long stats_bounds_avoided;
	static long stats_split_checks;
	static long stats_contain_checks;
	static long stats_split_full_checks;
	static long stats_split_checks_depths;
	ConvexPolygon():Polygon<D,T>(){

	}

	virtual ~ConvexPolygon(){};
	virtual ShapeType getType(){
		return CONVEX_POLYGON;
	}
	bool contains(const Point<D,T> & point);
	bool findContainingConvex(const Point<D,T> & point,NConvexPolygon<D,T> & polygon_out);

	bool containsInRange(const Point<D,T> & point, int firstVertex=0,int lastVertex=-1);
	bool containsInSplit(const Point<D,T> & point, int firstVertex=0,int lastVertex=-1);
	bool intersects(Shape<D,T> & s);


private:
	bool intersects2d(Shape<2,T> & s);
	bool edgesIntersectLine2d(LineSegment<2,T> & check, LineSegment<2,T> & intersection);
	bool containsInRange2d(const Point<2,T> & point, int firstVertex,int lastVertex);
	bool containsInSplit2d(const Point<2,T> & point, int firstVertex,int lastVertex, NConvexPolygon<2,T> & triangle_out);
	bool containsInSplit2d_helper(const Point<2,T> & point, int firstVertex,int lastVertex, NConvexPolygon<2,T> & triangle_out, int depth);
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
};

template<unsigned int D,class T>
long ConvexPolygon<D,T>::stats_triangle_avoided=0;
template<unsigned int D,class T>
long ConvexPolygon<D,T>::stats_bounds_avoided=0;
template<unsigned int D,class T>
long ConvexPolygon<D,T>::stats_split_checks=0;
template<unsigned int D,class T>
long ConvexPolygon<D,T>::stats_contain_checks=0;
template<unsigned int D,class T>
long ConvexPolygon<D,T>::stats_split_full_checks=0;

template<unsigned int D,class T>
long ConvexPolygon<D,T>::stats_split_checks_depths=0;



/**
 * A triangle defined by three vertices
 */
template<unsigned int D,class T>
class Triangle:public ConvexPolygon<D,T>{
public:
	//A line segment is defined by two end-points that it passes through
	Point<D,T> a;
	Point<D,T> b;
	Point<D,T> c;

	Triangle(){

	}
	Triangle(const Point<D,T> & a, const Point<D,T> & b,const Point<D,T> & c):a(a),b(b),c(c){

	}
	int size()const {
		return 3;
	}

	void update(){
		this->bounds_uptodate=false;
	}

	const Point<D,T>& operator [] (int index) const {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		if(index==0){
			return a;
		}else if (index==1){
			return b;
		}else{
			assert(index==2);
			return c;
		}
	}
	Point<D,T>&       operator [] (int index)       {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		if(index==0){
			return a;
		}else if (index==1){
			return b;
		}else{
			assert(index==2);
			return c;
		}
	}

};


template<unsigned int D,class T>
class NConvexPolygon:public ConvexPolygon<D,T>{
public:

	std::vector<Point<D,T>> vertices;


	NConvexPolygon():ConvexPolygon<D,T>(){

	}

	virtual ~NConvexPolygon(){};

	int size()const {
		return vertices.size();
	}

	void update(){
		if(!this->vertices_clockwise){
			reorderVertices();
		}
		this->bounds_uptodate=false;
	}

	void clear(){
		vertices.clear();
	}

	void addVertex(Point<D,T> p){
		this->vertices_clockwise=false;
		this->bounds_uptodate=false;
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
		this->bounds_uptodate=false;
	}

	void clearVertices(){
		vertices.clear();
		this->bounds_uptodate=false;
	}

	//Returns the vertices of the polygon, in clockwise order.

	std::vector<Point<D,T> > & getVertices(){
		if(!this->vertices_clockwise){
			reorderVertices();
		}
		this->dbg_orderClockwise();
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

	//put the vertices into clockwise order
	void reorderVertices(){
		if(D==2){
			reorderVertices2d();
		}else
			assert(false);
	}

private:
	// copy ops are private to prevent copying
	//Polygon(const Polygon& from); // no implementation
	NConvexPolygon& operator=(const NConvexPolygon& from); // no implementation

	//put the vertices into clockwise order
	void reorderVertices2d();

};







template<unsigned int D,class T>
bool ConvexPolygon<D,T>::findContainingConvex(const Point<D,T> & point,NConvexPolygon<D,T> & polygon_out){
	polygon_out.clear();
	stats_contain_checks++;
	if(this->size()==3){
		if(D==2){
			stats_triangle_avoided++;
			if(  this->pointInTriangle2d(point,(*this)[0],(*this)[1],(*this)[2])){
				polygon_out.addVertexUnchecked((*this)[0]);
				polygon_out.addVertexUnchecked((*this)[1]);
				polygon_out.addVertexUnchecked((*this)[2]);

				return true;
			}else{
				return false;
			}
		}
	}
	if(!Polygon<D,T>::boundContains(point)){
		assert(!containsInRange(point,0,this->size()-1));
		stats_bounds_avoided++;
		return false;
	}

	if(D==2){
		return containsInSplit2d((const Point<2,T> &) point, 0,this->size()-1,(NConvexPolygon<D,T> &)polygon_out);
	}else{
		assert(false);
	}
	return false;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::contains(const Point<D,T> & point)
{
	//stats_contain_checks++;
	static int iter = 0;
	if(++iter==309222){
		int a=1;
	}
	if( this->size()==3){
		if(D== 2)
			return this->pointInTriangle2d(point,(*this)[0],(*this)[1],(*this)[2]);
	}
	if(!this->boundContains(point)){
		//bool b = Polygon<D,T>::boundContains(point);
		assert(!containsInRange(point,0,this->size()-1));
		//stats_bounds_avoided++;
		return false;
	}

	if(Minisat::pipalg == PointInPolygonAlg::ALG_RECURSIVE_SPLIT){

		return containsInSplit(point,0,this->size()-1);
	}else{
		return containsInRange(point,0,this->size()-1);
	}
}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInRange(const Point<D,T> & point, int firstVertex,int lastVertex)
{
	if(D==2){
		return containsInRange2d((const Point<2,T> &) point, firstVertex,lastVertex);
	}else{
		assert(false);
	}
	return false;
}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInSplit(const Point<D,T> & point, int firstVertex,int lastVertex)
{
	if(D==2){
		static NConvexPolygon<2,T> ignore;
		return containsInSplit2d((const Point<2,T> &) point, firstVertex,lastVertex,ignore);
	}else{
		assert(false);
	}
	return false;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInSplit2d(const Point<2,T> & point, int firstVertex,int lastVertex, NConvexPolygon<2,T> & triangle_out){
	 stats_split_checks++;
	 //std::vector<Point<2,T> > &  w = (std::vector<Point<2,T> > & ) Polygon<D,T>::getVertices();
	 ConvexPolygon<2,T> & w = (ConvexPolygon<2,T> &) *this;
	 triangle_out.clear();
	 if(w.size()==0)
		 return false;
	 if(lastVertex<0)
		 lastVertex=w.size()-1;
	 assert(lastVertex>=0);
	 int n_verts = lastVertex-firstVertex+1;
	 if(lastVertex<firstVertex){
		 n_verts = lastVertex + ( w.size() - firstVertex );
	 }
	 assert(n_verts<=w.size());
	 this->dbg_orderClockwise();
	 if(n_verts==0)
		 return false;
	 else if(n_verts==1){
		 assert(lastVertex==firstVertex);
		 if( w[firstVertex]==point){
			 triangle_out.addVertex(w[firstVertex]);
			 return true;
		 }
		 return false;
	 }else if (n_verts==2){
		 assert(lastVertex!=firstVertex);
		 //from http://stackoverflow.com/a/11908158
		 //true if the point is between (inclusive) the other two points.
		 auto p1 = w[firstVertex];
		 auto p2 = w[lastVertex];
		 //check if the point lies on this line
		 if(crossDif(point, p1,p2)==0){
			 T dxl =p2.x-p1.x;
			 T dyl =p2.y-p1.y;
			 bool contains;
			 //check if the point is between the end points
			 if (abs(dxl) >= abs(dyl))
			    contains = dxl > 0 ?
					   p1.x <= point.x && point.x <= p2.x :
					   p2.x <= point.x && point.x <= p1.x;
			 else
			    contains = dyl > 0 ?
					   p1.y <= point.y && point.y <= p2.y :
					   p2.y <= point.y && point.y <= p1.y;

		     if(contains){
			    triangle_out.addVertex(w[firstVertex]);
			    triangle_out.addVertex(w[lastVertex]);
		     }
			 return contains;
		 }
		 return false;
	 }else if (n_verts==3){
		 bool contains= containsInRange(point,firstVertex,lastVertex);
		 if(contains){
			 triangle_out.addVertex(w[firstVertex]);
			 triangle_out.addVertex(w[firstVertex+1]);
			 triangle_out.addVertex(w[firstVertex+2]);
		 }
		 return contains;
	 }
	 stats_split_full_checks++;

	bool res= containsInSplit2d_helper(point, firstVertex,lastVertex, triangle_out,1);
	assert(res== containsInRange(point,firstVertex,lastVertex));
	return res;

}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInSplit2d_helper(const Point<2,T> & point,int first_vertex, int last_vertex, NConvexPolygon<2,T> & triangle_out, int depth){
	//recurse on this segment of the polygon, finding a triangle that contains the point.
	//precondition: the point is contained in this convex segment of the polygon

	//Noah's algorithm: pick 3 vertices in the polygon. 2 of them are adjacent, and the third is arbitrary (but should probably be the index that is farthest in both directions from the adjacent vertices)
	 //Check if they contain the point; if they do, return them.
	 //Else, check which of the two sides of the triangle with the non-adjacent vertex the point is on. Recurse on that sub polygon.

	 //When recursing, 2 of the three vertices are already selected (they are the vertices from the existing triangle), so we only have to pick one more vertex.
	 //Since we already know that the point isn't on the other side of those two vertices, we only have to check two sides in the case where the point is not contained.
	assert(first_vertex!=last_vertex);
/*
	static int iter = 0;
	if(++iter==66722){
		int a=1;
	}
	int it = iter;*/



	//std::vector<Point<2,T>> & polygon_vertices =(std::vector<Point<2,T>> & ) Polygon<D,T>::vertices;
	ConvexPolygon<2,T> & polygon_vertices = (ConvexPolygon<2,T> &)*this;
	Point<2,T> & a = polygon_vertices[first_vertex];
	Point<2,T> & b = polygon_vertices[last_vertex];
	assert(first_vertex<last_vertex);
	if(first_vertex+1==last_vertex) {
		stats_split_checks_depths+=depth;
		//Then this is a line. depending on our notion of containment, we either give up, or test if the line contains this point
		triangle_out.clear();
		triangle_out.addVertexUnchecked(a);
		triangle_out.addVertexUnchecked(b);
		if(triangle_out.contains(point)){
			return true;
		}else{
			triangle_out.clear();//give up
		}
		assert(!containsInRange(point,first_vertex,last_vertex));
		return false;
	}


	int mid_point = 0;
	if(first_vertex<last_vertex){
		mid_point = (last_vertex-first_vertex)/2 + first_vertex;
	}else{
		mid_point = (first_vertex-last_vertex)/2 + last_vertex;
	}
	Point<2,T> & c = polygon_vertices[mid_point];
	assert(mid_point != last_vertex);assert(mid_point!=first_vertex);
	/**/

	//compute the barycentric coordinates of this point
	auto & p0 = c;//b;
	auto & p1 = a;//c;
	auto & p2 = b;//a;


	assert(dbg_orderClockwise2dTri(p2,p1,p0));//intentionally reversing winding here, because the formula below is counter clockwise.
	T s = (p0.y*p2.x - p0.x*p2.y + (p2.y - p0.y)*point.x + (p0.x - p2.x)*point.y);
	T t = (p0.x*p1.y - p0.y*p1.x + (p0.y - p1.y)*point.x + (p1.x - p0.x)*point.y);
	T area2 =  (-p1.y*p2.x + p0.y*(-p1.x + p2.x) + p0.x*(p1.y - p2.y) + p1.x*p2.y); assert(area2>=0);
	/*T s = (a.y*b.x - a.x*b.y + (b.y - a.y)*point.x + (a.x - b.x)*point.y);
	T t = (a.x*c.y - a.y*c.x + (a.y - c.y)*point.x + (c.x - a.x)*point.y);*/
	//T gamma = 1-s-t;
	bool contained = ( s>=0 && t>=0 && (s+t<=area2));

	if(contained){
		triangle_out.clear();
		triangle_out.addVertexUnchecked(a);
		triangle_out.addVertexUnchecked(c);
		triangle_out.addVertexUnchecked(b);
		assert(triangle_out.contains(point));
		stats_split_checks_depths+=depth;
		return true;//we are done
	}else{
#ifndef NDEBUG
		ConvexPolygon<D,T> dbg_poly;
		dbg_poly.addVertexUnchecked(a);
		dbg_poly.addVertexUnchecked(c);
		dbg_poly.addVertexUnchecked(b);
		assert(!dbg_poly.contains(point));

#endif
/*

		Line<2,T> testLine(a,c);
		Line<2,T> testLine2(b,c);
		Line<2,T> testLine3(a,b);

		int as = testLine.whichSide(point);
		int ab = testLine2.whichSide(point);
		int ac = testLine3.whichSide(point);

		assert(testLine.whichSide(point)!= 0);//else we would have already found the point
*/
		if(t<0 && s>=0){
			//point is either between first_vertex,mid_point, or not in the triangle
			return containsInSplit2d_helper(point,first_vertex,mid_point,triangle_out,depth+1);
		}else if (s<0 && t>=0){
			//point is either between first_vertex,mid_point, or not in the triangle
			return containsInSplit2d_helper(point,mid_point,last_vertex,triangle_out,depth+1);
		}else{
			stats_split_checks_depths+=depth;
			//point is not contained.
			assert(!containsInRange(point,first_vertex,last_vertex));
			return false;//this point is not contained.

		}

	}
}


template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInRange2d(const Point<2,T> & point, int firstVertex,int lastVertex){
	//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
	//this is correct _only_ for convex polygons
	 //std::vector<Point<2,T> > &  w = (std::vector<Point<2,T> > & ) Polygon<D,T>::getVertices();
	 ConvexPolygon<2,T> & w = (ConvexPolygon<2,T> &) *this;
	 if(w.size()==0)
		 return false;
	 if(lastVertex<0)
		 lastVertex=w.size()-1;
	 assert(lastVertex>=0);
	 int n_verts = lastVertex-firstVertex+1;
	 if(lastVertex<firstVertex){
		 n_verts = lastVertex + ( w.size() - firstVertex );
	 }
	 assert(n_verts<=w.size());
	 this->dbg_orderClockwise();
	 if(n_verts==0)
		 return false;
	 else if(n_verts==1){
		 assert(lastVertex==firstVertex);
		 return w[firstVertex]==point;
	 }else if (n_verts==2){
		 assert(lastVertex!=firstVertex);
		 //from http://stackoverflow.com/a/11908158
		 //true if the point is between (inclusive) the other two points.
		 auto p1 = w[firstVertex];
		 auto p2 = w[lastVertex];
		 //check if the point lies on this line
		 if(crossDif(point, p1,p2)==0){
			 T dxl =p2.x-p1.x;
			 T dyl =p2.y-p1.y;

			 //check if the point is between the end points
			 if (abs(dxl) >= abs(dyl))
			   return dxl > 0 ?
					   p1.x <= point.x && point.x <= p2.x :
					   p2.x <= point.x && point.x <= p1.x;
			 else
			   return dyl > 0 ?
					   p1.y <= point.y && point.y <= p2.y :
					   p2.y <= point.y && point.y <= p1.y;
		 }
		 return false;
	 }

	 int endVertex = (lastVertex+1)% w.size();
		//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
		//this is correct _only_ for convex polygons
	 //note: this can also compute the area (which is the sum of p2[0]*p1[1] - p1[0]*p2[1]); could potentially combine these...
	 for(int n = 0;n<n_verts;n++){
		 int i = (firstVertex+n)% w.size();
		 Point<2,T> p1 = (i>0 ? w[i-1]:w.back()) - point;
		 Point<2,T> p2 = w[i]-point;
		 bool contained = (p2[0]*p1[1] - p1[0]*p2[1]) >=0;
		 if(!contained){
			 return false;
		 }
	 }
	 return true;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::intersects(Shape<D,T> & shape){
	if(D==2){
		return intersects2d((Shape<2,T> &) shape);
	}else{
		assert(false);
	}
	return false;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::edgesIntersectLine2d(LineSegment<2,T> & check,LineSegment<2,T> & out){
	// std::vector<Point<2,T> > &  w = this->getVertices();
	ConvexPolygon<2,T> & w = (ConvexPolygon<2,T>&)*this;
	for(int i = 0;i<w.size();i++){
		Point<2,T> & prev = i>0 ?  w[i-1]:w.back();
		Point<2,T> & p = w[i];
		out.a = prev;
		out.b = p;
		if(check.intersects(out)){
			return true;
		}
	}
	return false;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::intersects2d(Shape<2,T> & shape){
	if(shape.getType()==LINE_SEGMENT){
		LineSegment<2,T> & line = (LineSegment<2,T> &)shape;

		//first, check if either end point is contained
		if(this->contains2d(line.a) || this->contains2d(line.b))
			return true;
		static LineSegment<2,T>  ignore;
		//the line may still intersect even if neither end point is contained.
		//we could apply the SAT here. But instead, we're going to walk around the edges of the convex shape, and see if any of the edges intersect this line.
		return edgesIntersectLine2d(line,ignore);


	}else if(shape.getType()==CONVEX_POLYGON){
		ConvexPolygon<2,T> & c = (ConvexPolygon<2,T>&) shape;
		if(c.size()<this->size()){
			return c.intersects(*this);
		}

		if(this->size()==0 || c.size()==0){
			return false;
		}else if (this->size()==1){
			if(c.size()==1){
				//then the two vertices are considered to collide if they are identical
				return (*this)[0]==c[0];
			}else{
				return c.contains((*this)[0]);
			}
		}else if (c.size()==1){
			return contains(c[0]);
		}


		ConvexPolygon<2,T> &  w = (ConvexPolygon<2,T>&)*this;

		 //Separating Axis Theorem for collision detection between two convex polygons
		 //loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
		 //If any of the projections are non-intersection, then these don't collide; else, they do collide
		 for(int i = 0;i<this->size();i++){
			 auto & p = (*this)[i];
			 auto & prev = (*this)[i-1];
			 Point<2,T> edge = p-prev;
			 Point<2,T> un_normalized_normal(-edge.y, edge.x);

			 //now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
			 //Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
			 //(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
			 T left = std::numeric_limits<T>::infinity();
			 T right = -std::numeric_limits<T>::infinity();
			 for (auto & p:*this){
				 T projection = un_normalized_normal.dot(p);
				 if (projection < left) {
					  left = projection;
				 }
				 if (projection > right) {
					  right = projection;
				 }
			 }
			 bool overlaps = false;
			 bool seenLeft = false;
			 bool seenRight=true;
			 for (auto & p:c){
				 T projection = un_normalized_normal.dot(p);
				 if (projection >= left && projection <= right ) {
					 overlaps=true;
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
				 }
			 }
			 if(!overlaps && !(seenLeft&&seenRight)){
				 return false;
			 }
		 }

		 //now test the axis produced by the other polygon
		 for(int i = 0;i<c.size();i++){
			 auto & p = (*this)[i];
			 auto & prev = (*this)[i-1];
			 Point<2,T> edge = p-prev;
			 Point<2,T> un_normalized_normal(-edge.y, edge.x);

			 T left = std::numeric_limits<T>::infinity();
			 T right = -std::numeric_limits<T>::infinity();
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
			 bool seenRight=true;
			 bool overlaps = false;
			 for (auto & p:c){
				 T projection = un_normalized_normal.dot(p);
				 if (projection >= left && projection <= right ) {
					 overlaps=true;
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
				 }
			 }
			 if(!overlaps && !(seenLeft&&seenRight)){
				 return false;
			 }
		 }
		 //If no axis overlapped, then they did in fact intersect
		 return true;

	}
	assert(false);
	return false;
}


//put the vertices into clockwise order
template<unsigned int D,class T>
void NConvexPolygon<D,T>::reorderVertices2d(){
	this->vertices_clockwise=true;
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


#endif /* CONVEXPOLYGON_H_ */
