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
/**
 * A concrete, convex polygon (or, for D>2, a polytope)
 */
template<unsigned int D,class T>
class ConvexPolygon:public Polygon<D,T>{
public:

	virtual ~ConvexPolygon(){};
	virtual ShapeType getType(){
		return CONVEX_POLYGON;
	}
	bool contains(const Point<D,T> & point);
	bool containsInRange(const Point<D,T> & point, int firstVertex=0,int lastVertex=-1);
	bool containsInSplit(const Point<D,T> & point, int firstVertex=0,int lastVertex=-1);
	bool intersects(Shape<D,T> & s);


private:
	bool containsInRange2d(const Point<2,T> & point, int firstVertex,int lastVertex);
	bool containsInSplit2d(const Point<2,T> & point, int firstVertex,int lastVertex);
	bool containsInSplit2d_helper(const Point<2,T> & point, int firstVertex,int lastVertex, ConvexPolygon<2,T> & triangle_out);

};

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::contains(const Point<D,T> & point)
{
	static int iter = 0;
	if(++iter==309222){
		int a=1;
	}
	if(Polygon<D,T>::size()==3){
		if(D== 2)
			return Polygon<D,T>::pointInTriangle2d(point,Polygon<D,T>::vertices[0],Polygon<D,T>::vertices[1],Polygon<D,T>::vertices[2]);
	}
	if(!Polygon<D,T>::boundContains(point)){
		//bool b = Polygon<D,T>::boundContains(point);
		assert(!containsInRange(point,0,this->size()-1));
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
		return containsInSplit2d((const Point<2,T> &) point, firstVertex,lastVertex);
	}else{
		assert(false);
	}
	return false;
}

template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInSplit2d(const Point<2,T> & point, int firstVertex,int lastVertex){
	 std::vector<Point<2,T> > &  w = (std::vector<Point<2,T> > & ) Polygon<D,T>::getVertices();
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
	 Polygon<D,T>::dbg_orderClockwise();
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
	 }else if (n_verts==3){
		 return containsInRange(point,firstVertex,lastVertex);
	 }

	static ConvexPolygon<2,T> ignore;
	bool res= containsInSplit2d_helper(point, firstVertex,lastVertex, ignore);
	assert(res== containsInRange(point,firstVertex,lastVertex));
	return res;

}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInSplit2d_helper(const Point<2,T> & point,int first_vertex, int last_vertex, ConvexPolygon<2,T> & triangle_out){
	//recurse on this segment of the polygon, finding a triangle that contains the point.
	//precondition: the point is contained in this convex segment of the polygon

	//Noah's algorithm: pick 3 vertices in the polygon. 2 of them are adjacent, and the third is arbitrary (but should probably be the index that is farthest in both directions from the adjacent vertices)
	 //Check if they contain the point; if they do, return them.
	 //Else, check which of the two sides of the triangle with the non-adjacent vertex the point is on. Recurse on that sub polygon.

	 //When recursing, 2 of the three vertices are already selected (they are the vertices from the existing triangle), so we only have to pick one more vertex.
	 //Since we already know that the point isn't on the other side of those two vertices, we only have to check two sides in the case where the point is not contained.
	assert(first_vertex!=last_vertex);

	static int iter = 0;
	if(++iter==66722){
		int a=1;
	}
	int it = iter;



	std::vector<Point<2,T>> & polygon_vertices =(std::vector<Point<2,T>> & ) Polygon<D,T>::vertices;
	Point<2,T> & a = polygon_vertices[first_vertex];
	Point<2,T> & b = polygon_vertices[last_vertex];
	assert(first_vertex<last_vertex);
	if(first_vertex+1==last_vertex) {
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



	if(Polygon<D,T>::pointInTriangle2d(point,a,c,b)){
		triangle_out.clear();
		triangle_out.addVertexUnchecked(a);
		triangle_out.addVertexUnchecked(c);
		triangle_out.addVertexUnchecked(b);
		return true;//we are done
	}else{
		Line<2,T> testLine(a,c);
		assert(testLine.whichSide(point)!= 0);//else we would have already found the point

		if(testLine.whichSide(point)!= testLine.whichSide(b)){
			return containsInSplit2d_helper(point,first_vertex,mid_point,triangle_out);
		}else{
			testLine.a = b;
			assert(testLine.whichSide(point)!= 0);//else we would have already found the point
			if(testLine.whichSide(point)!= testLine.whichSide(a)){
				return containsInSplit2d_helper(point,mid_point,last_vertex,triangle_out);
			}else{
				assert(!containsInRange(point,first_vertex,last_vertex));
				return false;//this point is not contained.
			}


		}
	}
}


template<unsigned int D,class T>
bool ConvexPolygon<D,T>::containsInRange2d(const Point<2,T> & point, int firstVertex,int lastVertex){
	//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
	//this is correct _only_ for convex polygons
	 std::vector<Point<2,T> > &  w = (std::vector<Point<2,T> > & ) Polygon<D,T>::getVertices();
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
	 Polygon<D,T>::dbg_orderClockwise();
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



template<>
bool ConvexPolygon<2,mpq_class>::intersects(Shape<2,mpq_class> & s);

template<>
bool ConvexPolygon<2,double>::intersects(Shape<2,double> & s);


template<>
bool ConvexPolygon<1,double>::intersects(Shape<1,double> & s);


#endif /* CONVEXPOLYGON_H_ */
