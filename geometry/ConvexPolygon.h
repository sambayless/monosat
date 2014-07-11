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
#include <iostream>
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
	bool contains(const Point<D,T> & point, bool inclusive);
	bool contains(const Point<D,T> & point,NConvexPolygon<D,T> & polygon_out, bool inclusive);

	bool containsInRange(const Point<D,T> & point, int firstVertex,int lastVertex, bool inclusive);
	bool containsInSplit(const Point<D,T> & point, int firstVertex,int lastVertex, bool inclusive);
	bool intersects(Shape<D,T> & s, bool inclusive);

	bool intersects(Shape<D,T> & s, NConvexPolygon<D,T> & polygon_out, bool inclusive);

private:
	bool intersects2d(Shape<2,T> & s, bool inclusive);
	bool intersects2d(Shape<2,T> & s, NConvexPolygon<2,T> & polygon_out, bool inclusive);
	bool edgesIntersectLine2d(LineSegment<2,T> & check, LineSegment<2,T> & intersection, bool inclusive);
	bool containsInRange2d(const Point<2,T> & point, int firstVertex,int lastVertex, bool inclusive);
	bool containsInSplit2d(const Point<2,T> & point, int firstVertex,int lastVertex, NConvexPolygon<2,T> & polygon_out, bool inclusive);
	bool containsInSplit2d_helper(const Point<2,T> & point, int firstVertex,int lastVertex, NConvexPolygon<2,T> & polygon_out, int depth, bool inclusive);
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
		assert(this->dbg_orderClockwise());
		assert(this->dbg_boundsUpToDate());
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
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		return vertices[index];
	}
	Point<D,T>&       operator [] (int index)       {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
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
std::ostream & operator<<(std::ostream & str, ConvexPolygon<D,T>  & polygon){
	str << "ConvexPolygon=[";
	for (const auto & p:polygon){
		str<<p <<",";
	}
	str<<"]";
	return str;
}

template<unsigned int D,class T>
std::ostream & operator<<(std::ostream & str, Triangle<D,T>  & polygon){
	str << "Triangle=[";
	for (const auto & p:polygon){
		str<<p <<",";
	}
	str<<"]";
	return str;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::contains(const Point<D,T> & point,NConvexPolygon<D,T> & polygon_out, bool inclusive){
	polygon_out.clear();
	stats_contain_checks++;
	if(this->size()==3){
		if(D==2){
			stats_triangle_avoided++;
			if(  this->pointInTriangle2d(point,(*this)[0],(*this)[1],(*this)[2],inclusive)){
				polygon_out.addVertexUnchecked((*this)[0]);
				polygon_out.addVertexUnchecked((*this)[1]);
				polygon_out.addVertexUnchecked((*this)[2]);

				return true;
			}else{
				return false;
			}
		}
	}
	if(!this->boundContains(point,inclusive)){
		assert(!containsInRange(point,0,this->size()-1,inclusive));
		stats_bounds_avoided++;
		return false;
	}

	if(D==2){
		return containsInSplit2d((const Point<2,T> &) point, 0,this->size()-1,(NConvexPolygon<D,T> &)polygon_out,inclusive);
	}else{
		assert(false);
	}
	return false;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::contains(const Point<D,T> & point, bool inclusive)
{
	//stats_contain_checks++;
	static int iter = 0;
	if(++iter==309222){
		int a=1;
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
	}
}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInRange(const Point<D,T> & point, int firstVertex,int lastVertex, bool inclusive)
{
	this->update();
	if(D==2){
		return containsInRange2d((const Point<2,T> &) point, firstVertex,lastVertex,inclusive);
	}else{
		assert(false);
	}
	return false;
}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInSplit(const Point<D,T> & point,  int firstVertex,int lastVertex,bool inclusive)
{
	this->update();
	if(D==2){
		static NConvexPolygon<2,T> ignore;
		return containsInSplit2d((const Point<2,T> &) point,firstVertex,lastVertex,ignore,inclusive);
	}else{
		assert(false);
	}
	return false;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInSplit2d(const Point<2,T> & point, int firstVertex,int lastVertex, NConvexPolygon<2,T> & polygon_out, bool inclusive){
	 stats_split_checks++;
	 //std::vector<Point<2,T> > &  w = (std::vector<Point<2,T> > & ) Polygon<D,T>::getVertices();
	 ConvexPolygon<2,T> & w = (ConvexPolygon<2,T> &) *this;
	 polygon_out.clear();
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
		 if(!inclusive)
			 return false;
		 assert(lastVertex==firstVertex);
		 if( w[firstVertex]==point){
			 polygon_out.addVertex(w[firstVertex]);
			 return true;
		 }
		 return false;
	 }else if (n_verts==2){
		 if(!inclusive)
			 return false;
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
		    	 polygon_out.addVertex(w[firstVertex]);
		    	 polygon_out.addVertex(w[lastVertex]);
		     }
			 return contains;
		 }
		 return false;
	 }else if (n_verts==3){
		 bool contains= containsInRange(point,firstVertex,lastVertex, inclusive);
		 if(contains){
			 polygon_out.addVertex(w[firstVertex]);
			 polygon_out.addVertex(w[firstVertex+1]);
			 polygon_out.addVertex(w[firstVertex+2]);
		 }
		 return contains;
	 }
	 stats_split_full_checks++;

	bool res= containsInSplit2d_helper(point, firstVertex,lastVertex, polygon_out,1,inclusive);
	if(res && !inclusive){
		//do an extra check to correct for edge cases.
		assert(polygon_out.contains(point,true));
		if(!polygon_out.contains(point,false)){
			//then there are two possibilities. Either the point is exactly on an edge of the polygon (in whch case it is not contained unless we are being inclsuvie)
			//OR, we can expan the triangle into a quadrilateral that _does_ contain the point.
			//so, check which of the three edges of the triangle the point is lying on (if it is on a vertex of the triangle, then it must be of the border of the polygon, so it is not included)
			assert(polygon_out.size()==3);
			int containing_edge=-1;
			for(int i = 0;i<polygon_out.size();i++){
				Point<2,T> & prev = polygon_out[i-1];
				Point<2,T> & p = polygon_out[i];
				if(eq_epsilon(crossDif(prev,p,point))){
					containing_edge=i;
					break;
				}
			}
			assert(containing_edge>=0);
			int pIndex = -1;
			int prevIndex = -1;
			int nIndex = -1;
			//now check whether that edge is on the border of the polygon
			for(int i = 0;i<this->size();i++){
				Point<2,T> & prev = (*this)[i-1];
				Point<2,T> & p = (*this)[i];
				if(p.getID()== polygon_out[containing_edge].getID()){
						pIndex=i;
				}else if (p.getID()==polygon_out[containing_edge-1].getID()){
					prevIndex=i;
				}else if (p.getID()==polygon_out[containing_edge+1].getID()){
					nIndex=i;
				}
			}
			if(pIndex<prevIndex){
				std::swap(pIndex,prevIndex);
			}
			if((pIndex-prevIndex) % this->size()==1){
				//The containing edge is on the border of the polygon, so the point is NOT contained exclusively within the polygon.
				res=false;
				polygon_out.clear();

			}else{
				//None of the edges of the hull contain the point, so the containing edge is interior to the polygon - meaning that the point is exclusively contained in the polygon, but that the current containing triangle
				//either needs to be swapped for some other triangle, or widened to a quadrilateral.
				//(since there exist edge cases where NO internal triangle will exclusively contain the point, we will just widen to a quadrilateral, here.)
				assert(prevIndex<pIndex);
				assert(pIndex>-1);
				assert(prevIndex>-1);
				assert(nIndex>-1);
				assert(pIndex!=prevIndex);
				assert(pIndex!=nIndex);

				if(prevIndex<nIndex && nIndex<pIndex){
					//ok, pick any vertex outside of pIndex, prevIndex in order to expand the quadrilateral to contain point
					assert((nIndex+1) % this->size() != prevIndex);
					polygon_out.addVertex((*this)[pIndex+1]);
				}else{
					assert(nIndex<prevIndex || nIndex>pIndex);
					assert((pIndex-1) % this->size() != prevIndex);
					polygon_out.addVertex((*this)[pIndex-1]);
				}
				assert(polygon_out.containsInRange(point,0,polygon_out.size()-1,true));
				assert(polygon_out.containsInRange(point,0,polygon_out.size()-1,false));
			}

		}
	}

	assert(res== containsInRange(point,firstVertex,lastVertex,inclusive));
	return res;

}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInSplit2d_helper(const Point<2,T> & point,int first_vertex, int last_vertex, NConvexPolygon<2,T> & triangle_out, int depth, bool inclusive){
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
		if(!inclusive){
			return false;
		}
		//Then this is a line. depending on our notion of containment, we either give up, or test if the line contains this point
		triangle_out.clear();
		triangle_out.addVertexUnchecked(a);
		triangle_out.addVertexUnchecked(b);
		if(triangle_out.contains(point,inclusive)){
			return true;
		}else{
			triangle_out.clear();//give up
		}
		assert(!containsInRange(point,first_vertex,last_vertex, inclusive));
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
	bool contained;

	//checks for inclusiveness will be done outside this method
	//if(inclusive)
	contained= ( s>=0 && t>=0 && (s+t<=area2));
	//else
	//	contained= ( s>0 && t>0 && (s+t<area2));

	if(contained){
		triangle_out.clear();
		triangle_out.addVertexUnchecked(a);
		triangle_out.addVertexUnchecked(c);
		triangle_out.addVertexUnchecked(b);
		assert(triangle_out.contains(point,true));
		stats_split_checks_depths+=depth;
		return true;//we are done
	}else{
#ifndef NDEBUG
		NConvexPolygon<D,T> dbg_poly;
		dbg_poly.addVertexUnchecked(a);
		dbg_poly.addVertexUnchecked(c);
		dbg_poly.addVertexUnchecked(b);
		assert(!dbg_poly.contains(point,true));

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
			return containsInSplit2d_helper(point,first_vertex,mid_point,triangle_out,depth+1,inclusive);
		}else if (s<0 && t>=0){
			//point is either between first_vertex,mid_point, or not in the triangle
			return containsInSplit2d_helper(point,mid_point,last_vertex,triangle_out,depth+1,inclusive);
		}else{
			stats_split_checks_depths+=depth;
			//point is not contained.
			assert(!containsInRange(point,first_vertex,last_vertex,inclusive));
			return false;//this point is not contained.

		}

	}
}


template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInRange2d(const Point<2,T> & point, int firstVertex,int lastVertex, bool inclusive){
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
		 return inclusive && w[firstVertex]==point;
	 }else if (n_verts==2){
		 if (!inclusive)
			 return false;
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
		 bool contained;
		 if(inclusive)
			 contained= (p2[0]*p1[1] - p1[0]*p2[1]) >=0;
		 else
			 contained= (p2[0]*p1[1] - p1[0]*p2[1]) >0;
		 if(!contained){
			 return false;
		 }
	 }
	 return true;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::intersects(Shape<D,T> & shape,NConvexPolygon<D,T> & out, bool inclusive){
	out.clear();
	if(D==2){
		return intersects2d((Shape<2,T> &) shape,(NConvexPolygon<2,T> &)out,inclusive);
	}else{
		assert(false);
	}
	return false;
}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::intersects(Shape<D,T> & shape, bool inclusive){
	if(D==2){
		return intersects2d((Shape<2,T> &) shape,inclusive);
	}else{
		assert(false);
	}
	return false;
}


template<unsigned int D,class T>
bool ConvexPolygon<D,T>::edgesIntersectLine2d(LineSegment<2,T> & check,LineSegment<2,T> & out, bool inclusive){
	// std::vector<Point<2,T> > &  w = this->getVertices();
	ConvexPolygon<2,T> & w = (ConvexPolygon<2,T>&)*this;
	for(int i = 0;i<w.size();i++){
		Point<2,T> & prev = i>0 ?  w[i-1]:w.back();
		Point<2,T> & p = w[i];
		out.a = prev;
		out.b = p;
		if(check.intersects(out,inclusive)){
			return true;
		}
	}
	return false;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::intersects2d(Shape<2,T> & shape, bool inclusive){
	if(shape.getType()==LINE_SEGMENT){
		LineSegment<2,T> & line = (LineSegment<2,T> &)shape;

		//first, check if either end point is contained
		if(this->contains(line.a,inclusive) || this->contains(line.b,inclusive))
			return true;
		static LineSegment<2,T>  ignore;
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
			 bool overlaps = false;
			 bool seenLeft = false;
			 bool seenRight=false;
			 for (auto & p:c){
				 T projection = un_normalized_normal.dot(p);
				 if(inclusive){
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
					 }else if (seenLeft && projection > left ){
						 seenRight=true;
						 break;
					 }else if (seenRight && projection < right ){
						 seenRight=true;
						 break;
					 }
				 }else{
					 if (projection > left && projection < right ) {
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
					 }else if (seenLeft && projection > left ){
						 seenRight=true;
						 break;
					 }else if (seenRight && projection < right ){
						 seenRight=true;
						 break;
					 }
				 }
			 }
			 if(!overlaps && !(seenLeft&&seenRight)){
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
			 bool overlaps = false;
			 for (auto & p:*this){
				 T projection = un_normalized_normal.dot(p);
				 if(inclusive){
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
					 }else if (seenLeft && projection > left ){
						 seenRight=true;
						 break;
					 }else if (seenRight && projection < right ){
						 seenRight=true;
						 break;
					 }
				 }else{
					 if (projection > left && projection < right ) {
						 overlaps=true;
						 break;
					 }else if (projection <= left ){
						 seenLeft=true;
						 if(seenRight){
							 break;
						 }
					 }else if (projection >= right){
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
				 }
			 }
			 if(!overlaps && !(seenLeft&&seenRight)){
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

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::intersects2d(Shape<2,T> & shape, NConvexPolygon<2,T> & out, bool inclusive){
	if(this->size()==0)
		return false;
	if(shape.getType()==LINE_SEGMENT){
		LineSegment<2,T> & line = (LineSegment<2,T> &)shape;

		//first, check if either end point is contained
		if(this->contains(line.a,out,inclusive) || this->contains(line.b,out,inclusive))
			return true;
		if(this->size()==1){
			return line.contains((*this)[0],inclusive);
		}
		static LineSegment<2,T>  store;
		//the line may still intersect even if neither end point is contained.
		//we could apply the SAT here. But instead, we're going to walk around the edges of the convex shape, and see if any of the edges intersect this line.
		bool r = edgesIntersectLine2d(line,store,inclusive);
		if(r){
			out.addVertex(store.a);
			out.addVertex(store.b);

		}
#ifndef NDEBUG
		NConvexPolygon<2,T> test;
		test.addVertex(line.a);
		test.addVertex(line.b);
		if(this->intersects2d(test,inclusive)!=r){
			std::cout<<line<<"\n";
			std::cout<<(*this)<<"\n";
			edgesIntersectLine2d(line,store,inclusive);
		}
		assert(this->intersects2d(test,inclusive)==r);
#endif
		return r;


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

	assert(this->dbg_orderClockwise());
}


#endif /* CONVEXPOLYGON_H_ */
