/*
 * Polygon.cpp
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#include "ConvexPolygon.h"
#include "mtl/Vec.h"
template<>
bool ConvexPolygon<2,mpq_class>::containsInRange(const Point<2,mpq_class> & point, int firstVertex,int lastVertex){
	//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
	//this is correct _only_ for convex polygons
	 std::vector<Point<2,mpq_class> > &  w = getVertices();
	 dbg_orderClockwise();
	 if(w.size()==0)
		 return false;
	 else if(w.size()==1){
		 return w[0]==point;
	 }/*else if (w.size()==2){
		 //true if the point is between (inclusive) the other two points.

	 }*/
	 if(lastVertex<0)
		 lastVertex=w.size()-1;
	 //note: this can also compute the area (which is the sum of p2[0]*p1[1] - p1[0]*p2[1]); could potentially combine these...
	 for(int i = firstVertex;i!=lastVertex;i=(i+1% w.size() )){
		 Point<2,mpq_class> p1 = i>0 ? (w[i-1]-point):((w[lastVertex]-point));
		 Point<2,mpq_class> p2 = w[i]-point;
		 mpq_class t =  (p2[0]*p1[1] - p1[0]*p2[1]);
		 bool contained = (p2[0]*p1[1] - p1[0]*p2[1]) >=0;
		 if(!contained){
			 return false;
		 }
	 }
	 return true;
}

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
			 double left = std::numeric_limits<double>::max();
			 double right = std::numeric_limits<double>::min();
			 for (auto & p:*this){
				 double projection = un_normalized_normal.dot(p);
				 if (projection < left) {
					  left = projection;
				 } else if (projection > right) {
					  right = projection;
				 }
			 }
			 bool overlaps = false;
			 for (auto & p:c){
				 double projection = un_normalized_normal.dot(p);
				 if (projection >= left && projection <= right ) {
					 overlaps=true;
					 break;
				 }
			 }
			 if(!overlaps){
				 return false;
			 }
		 }

		 //now test the axis produced by the other polygon
		 for(int i = 0;i<c.size();i++){
			 auto & p = (*this)[i];
			 auto & prev = (*this)[i-1];
			 Point<2,double> edge = p-prev;
			 Point<2,double> un_normalized_normal(-edge.y, edge.x);

			 double left = std::numeric_limits<double>::max();
			 double right = std::numeric_limits<double>::min();
			 for (auto & p:*this){
				 double projection = un_normalized_normal.dot(p);
				 if (projection < left) {
					  left = projection;
				 } else if (projection > right) {
					  right = projection;
				 }
			 }
			 bool overlaps = false;
			 for (auto & p:c){
				 double projection = un_normalized_normal.dot(p);
				 if (projection >= left && projection <= right ) {
					 overlaps=true;
					 break;
				 }
			 }
			 if(!overlaps){
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
		ConvexPolygon<2,mpq_class> & c = *(ConvexPolygon<2,mpq_class>*) &shape;
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
			 mpq_class left = std::numeric_limits<mpq_class>::max();
			 mpq_class right = std::numeric_limits<mpq_class>::min();
			 for (auto & p:*this){
				 mpq_class projection = un_normalized_normal.dot(p);
				 if (projection < left) {
					  left = projection;
				 } else if (projection > right) {
					  right = projection;
				 }
			 }
			 bool overlaps = false;
			 for (auto & p:c){
				 mpq_class projection = un_normalized_normal.dot(p);
				 if (projection >= left && projection <= right ) {
					 overlaps=true;
					 break;
				 }
			 }
			 if(!overlaps){
				 return false;
			 }
		 }

		 //now test the axis produced by the other polygon
		 for(int i = 0;i<c.size();i++){
			 auto & p = (*this)[i];
			 auto & prev = (*this)[i-1];
			 Point<2,mpq_class> edge = p-prev;
			 Point<2,mpq_class> un_normalized_normal(-edge.y, edge.x);

			 mpq_class left = std::numeric_limits<mpq_class>::max();
			 mpq_class right = std::numeric_limits<mpq_class>::min();
			 for (auto & p:*this){
				 mpq_class projection = un_normalized_normal.dot(p);
				 if (projection < left) {
					  left = projection;
				 } else if (projection > right) {
					  right = projection;
				 }
			 }
			 bool overlaps = false;
			 for (auto & p:c){
				 mpq_class projection = un_normalized_normal.dot(p);
				 if (projection >= left && projection <= right ) {
					 overlaps=true;
					 break;
				 }
			 }
			 if(!overlaps){
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
bool ConvexPolygon<2,mpq_class>::contains(const Point<2,mpq_class> & point){
	return containsInRange(point,0,this->size()-1);
}
template<>
bool ConvexPolygon<2,double>::containsInRange(const Point<2,double> & point, int firstVertex,int lastVertex){
	//From http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
	//this is correct _only_ for convex polygons
	 std::vector<Point<2,double> > &  w = getVertices();
	 dbg_orderClockwise();
	 if(w.size()==0)
		 return false;
	 else if(w.size()==1){
		 return w[0]==point;
	 }/*else if (w.size()==2){
		 //true if the point is between (inclusive) the other two points.

	 }*/
	 if(lastVertex<0)
		 lastVertex=w.size()-1;
	 //note: this can also compute the area (which is the sum of p2[0]*p1[1] - p1[0]*p2[1]); could potentially combine these...
	 for(int i = firstVertex;i!=lastVertex;i=(i+1% w.size() )){
		 Point<2,double> p1 = i>0 ? (w[i-1]-point):((w[lastVertex]-point));
		 Point<2,double> p2 = w[i]-point;
		 bool contained = (p2[0]*p1[1] - p1[0]*p2[1]) >=0;
		 if(!contained){
			 return false;
		 }
	 }
	 return true;
}



template<>
bool ConvexPolygon<1,double>::containsInRange(const Point<1,double> & point, int firstVertex,int lastVertex){
	return false;
}
template<>
bool ConvexPolygon<1,double>::contains(const Point<1,double> & point){
	return containsInRange(point,0,this->size()-1);
}
template<>
bool ConvexPolygon<2,double>::contains(const Point<2,double> & point){
	return containsInRange(point,0,this->size()-1);
}
template<unsigned int D,class T>
bool ConvexPolygon<D,T>::contains(const Point<D,T> & point)
{
	return containsInRange(point,0,this->size()-1);
}


template<>
bool ConvexPolygon<1,double>::intersects(Shape<1,double> & shape){
	return false;
}
