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

/**
 * A concrete polygon (or, for D>2, a polytope)
 */
template<unsigned int D,class T=double>
class Polygon:public Shape<D>{
public:
	//List of vertices in clockwise order
	std::vector<Point<D,T>> vertices;
	Point<D,T> circleCenter;//should generalize this to an arbitrary bounding volume...
	T circleRadius;
	bool vertices_clockwise=false;
	virtual ~Polygon(){};
	virtual ShapeType getType(){
		return POLYGON;
	}
	virtual bool contains(const Point<D,T> & point);

	virtual bool intersects(Shape<D,T> & s){
		return false;
	}

	int size()const {
		return vertices.size();
	}

	void update(){
		if(!vertices_clockwise){
			reorderVertices();
		}
		updateCircleBound();
	}

	void clear(){
		vertices.clear();
	}

	int size(){
		return vertices.size();
	}

	void addVertex(Point<D,T> & p){
		vertices_clockwise=false;
		vertices.push_back(p);
	}
	//add a vertex, assuming that it will preserve clockwise order
	void addVertexUnchecked(Point<D,T> & p){
		vertices.push_back(p);
		assert(dbg_orderClockwise());
	}

	void popVertex(){
		vertices.pop();
	}

	void clearVertices(){
		vertices.clear();
	}

	//Returns the vertices of the polygon, in clockwise order.

	std::vector<Point<D,T> > & getVertices(){
		if(!vertices_clockwise){
			reorderVertices();
		}
		return vertices;
	}


	void updateCircleBound();
	virtual T getArea();
	virtual T getPerimeter();
	private:

	bool dbg_orderClockwise(){
		return true;
	}

	void reorderVertices();



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


	    wrap_iterator( std::vector<ValueType> & verts, int pos):verts(verts){}; // private constructor for begin, end
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
bool Polygon<2,double>::contains(const Point<2,double> & point);

template<>
void Polygon<2,double>::reorderVertices();

template<>
double Polygon<2,double>::getArea();

template<>
double Polygon<2,double>::getPerimeter();

template<unsigned int D,class T>
void Polygon<D,T>::updateCircleBound(){
	circleCenter.zero();
	for(int i = 0;i<vertices.size();i++){
		circleCenter+=vertices[i];
	}
	circleCenter/=T(vertices.size());
	circleRadius=T(0);
	for(int i = 0;i<vertices.size();i++){
		T dist = circleCenter.distance( vertices[i]);
		if(dist>circleRadius){
			circleRadius=dist;
		}
	}
}

#endif /* POLYGON_H_ */
