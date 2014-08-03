/*
 * CGALConvexPolygon.h
 *
 *  Created on: Aug 2, 2014
 *      Author: sam
 */

#ifndef CGALCONVEXPOLYGON_H_
#define CGALCONVEXPOLYGON_H_

#include "../ConvexPolygon.h"

#include <CGAL/basic.h>
// GMP is installed. Use the GMP rational number-type.
#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
typedef CGAL::Cartesian<CGAL::Gmpq>                Kernel;
typedef Kernel::Point_2                             Point_2;
typedef CGAL::Polygon_2<Kernel>                     Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>          Polygon_with_holes_2;

template<class T>
class CGALConvexPolygon2:public ConvexPolygon<2,T>{
public:

	Polygon_2 polygon;

	CGALConvexPolygon2():ConvexPolygon<2,T>(){

	}

	virtual ~CGALConvexPolygon2(){};

	int size()const {
		return polygon.size();
	}

	void update(){

		this->bounds_uptodate=false;
	}

	void clear(){
		polygon.clear();
		this->bounds_uptodate=false;
	}

	void addVertex(Point<2,T> p){
		polygon.push_back(Point_2(p[0],p[1]));
		this->bounds_uptodate=false;

	}

	//add a vertex, assuming that it will preserve clockwise order
	void addVertexUnchecked(Point<2,T>  p){
		polygon.push_back(Point_2(p[0],p[1]));
	}

	void popVertex(){
		polygon.erase(--polygon.vertices_end());
		this->bounds_uptodate=false;
	}

	//Returns the vertices of the polygon, in clockwise order.

/*	std::vector<Point<2,T> > & getVertices(){
		if(!this->vertices_clockwise){
			reorderVertices();
		}
		this->dbg_orderClockwise();
		return vertices;
	}*/
	const Point<2,T>& operator [] (int index) const {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		return Point<2,T>(polygon[index]);
		//return vertices[index];
	}
	Point<2,T>&       operator [] (int index)       {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		return Point<2,T>(polygon[index]);
	}



private:
	// copy ops are private to prevent copying
	//Polygon(const Polygon& from); // no implementation
	CGALConvexPolygon2& operator=(const CGALConvexPolygon2& from); // no implementation


};
/*
template<class T>
class CGALConvexPolygonWithHoles2:public ConvexPolygon<2,T>{
public:

	Polygon_with_holes_2 polygon;

	CGALConvexPolygonWithHoles2():ConvexPolygon<2,T>(){

	}

	virtual ~CGALConvexPolygonWithHoles2(){};

	int size()const {
		return polygon.size();
	}

	void update(){

		this->bounds_uptodate=false;
	}

	void clear(){
		polygon.clear();
		this->bounds_uptodate=false;
	}

	void addVertex(Point<2,T> p){
		polygon.push_back(Point_2(p[0],p[1]));
		this->bounds_uptodate=false;

	}

	//add a vertex, assuming that it will preserve clockwise order
	void addVertexUnchecked(Point<2,T>  p){
		polygon.push_back(Point_2(p[0],p[1]));
	}

	void popVertex(){
		polygon.erase(--polygon.vertices_end());
		this->bounds_uptodate=false;
	}

	//Returns the vertices of the polygon, in clockwise order.

	std::vector<Point<2,T> > & getVertices(){
		if(!this->vertices_clockwise){
			reorderVertices();
		}
		this->dbg_orderClockwise();
		return vertices;
	}
	const Point<2,T>& operator [] (int index) const {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		return Point<2,T>(polygon[index]);
		//return vertices[index];
	}
	Point<2,T>&       operator [] (int index)       {
		index = index %size();
		if(index<0){
			index+=size();
		}
		assert(index>=0);assert(index<size());
		return Point<2,T>(polygon[index]);
	}



private:
	// copy ops are private to prevent copying
	//Polygon(const Polygon& from); // no implementation
	CGALConvexPolygonWithHoles2& operator=(const CGALConvexPolygonWithHoles2& from); // no implementation


};*/

#endif /* CGALCONVEXPOLYGON_H_ */
