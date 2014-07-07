/*
 * Polygon.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#include "ConvexPolygon.h"
#include "mtl/Vec.h"
#include <cmath>


template<>
bool ConvexPolygon<2,double>::intersects(Shape<2,double> & shape){
	if(shape.getType()==CONVEX_POLYGON){
		ConvexPolygon<2,double> & c = *(ConvexPolygon<2,double> *) &shape;
		if(c.size()<size()){
			return c.intersects(*this);
		}

		if(size()==0 || c.size()==0){
			return false;
		}else if (size()==1){
			if(c.size()==1){
				//then the two vertices are considered to collide if they are identical
				return (*this)[0]==c[0];
			}else{
				return c.contains((*this)[0]);
			}
		}else if (c.size()==1){
			return contains(c[0]);
		}


		 std::vector<Point<2,double> > &  w = getVertices();

		 //Separating Axis Theorem for collision detection between two convex polygons
		 //loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
		 //If any of the projections are non-intersection, then these don't collide; else, they do collide
		 for(int i = 0;i<size();i++){
			 auto & p = (*this)[i];
			 auto & prev = (*this)[i-1];
			 Point<2,double> edge = p-prev;
			 Point<2,double> un_normalized_normal(-edge.y, edge.x);

			 //now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
			 //Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
			 //(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
			 double left = std::numeric_limits<double>::infinity();
			 double right = -std::numeric_limits<double>::infinity();
			 for (auto & p:*this){
				 double projection = un_normalized_normal.dot(p);
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
				 double projection = un_normalized_normal.dot(p);
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
			 Point<2,double> edge = p-prev;
			 Point<2,double> un_normalized_normal(-edge.y, edge.x);

			 double left = std::numeric_limits<double>::infinity();
			 double right = -std::numeric_limits<double>::infinity();
			 for (auto & p:*this){
				 double projection = un_normalized_normal.dot(p);
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
				 double projection = un_normalized_normal.dot(p);
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

template<>
bool ConvexPolygon<2,mpq_class>::intersects(Shape<2,mpq_class> & shape){
	if(shape.getType()==CONVEX_POLYGON){
		ConvexPolygon<2,mpq_class> & c = *(ConvexPolygon<2,mpq_class> *) &shape;
		if(c.size()<size()){
			return c.intersects(*this);
		}

		if(size()==0 || c.size()==0){
			return false;
		}else if (size()==1){
			if(c.size()==1){
				//then the two vertices are considered to collide if they are identical
				return (*this)[0]==c[0];
			}else{
				return c.contains((*this)[0]);
			}
		}else if (c.size()==1){
			return contains(c[0]);
		}


		 std::vector<Point<2,mpq_class> > &  w = getVertices();

		 //Separating Axis Theorem for collision detection between two convex polygons
		 //loop through each edge in _each_ polygon and project both polygons onto that edge's normal.
		 //If any of the projections are non-intersection, then these don't collide; else, they do collide
		 for(int i = 0;i<size();i++){
			 auto & p = (*this)[i];
			 auto & prev = (*this)[i-1];
			 Point<2,mpq_class> edge = p-prev;
			 Point<2,mpq_class> un_normalized_normal(-edge.y, edge.x);

			 //now project both polygons onto to this normal and see if they overlap, by finding the minimum and maximum distances
			 //Note that since we are NOT normalizing the normal vector, the projection is distorted along that vector
			 //(this still allows us to check overlaps, but means that the minimum distance found between the two shapes may be incorrect)
			 mpq_class left = std::numeric_limits<mpq_class>::infinity();
			 mpq_class right = -std::numeric_limits<mpq_class>::infinity();
			 for (auto & p:*this){
				 mpq_class projection = un_normalized_normal.dot(p);
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
				 mpq_class projection = un_normalized_normal.dot(p);
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
			 Point<2,mpq_class> edge = p-prev;
			 Point<2,mpq_class> un_normalized_normal(-edge.y, edge.x);

			 mpq_class left = std::numeric_limits<mpq_class>::infinity();
			 mpq_class right = -std::numeric_limits<mpq_class>::infinity();
			 for (auto & p:*this){
				 mpq_class projection = un_normalized_normal.dot(p);
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
				 mpq_class projection = un_normalized_normal.dot(p);
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

template<>
bool ConvexPolygon<2,double>::containsInRange(const Point<2,double> & point, int firstVertex,int lastVertex){

	 std::vector<Point<2,double> > &  w = getVertices();
	 if(w.size()==0)
		 return false;
	 if(lastVertex<0)
		 lastVertex+=w.size();
	 assert(lastVertex>=0);
	 int n_verts = lastVertex-firstVertex+1;
	 if(lastVertex<firstVertex){
		 n_verts = lastVertex + ( w.size() - firstVertex );
	 }
	 assert(n_verts<=w.size());
	 dbg_orderClockwise();
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
			 double dxl =p2.x-p1.x;
			 double dyl =p2.y-p1.y;

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
		 Point<2,double> p1 = i>0 ? (w[i-1]-point):((w.back()-point));
		 Point<2,double> p2 = w[i]-point;
		 bool contained = (p2[0]*p1[1] - p1[0]*p2[1]) >=0;
		 if(!contained){
			 return false;
		 }
	 }
	 return true;
}





template<>
bool ConvexPolygon<1,double>::intersects(Shape<1,double> & shape){
	return false;
}
