
#ifndef DELAUNAYPOLYPARTITION_H_
#define DELAUNAYPOLYPARTITION_H_


#include "Delaunay.h"

#include "PolygonSet.h"
#include <gmpxx.h>
#include <vector>

//interface to the polypartition library
template<unsigned int D,class T>
class DelaunayPolypartition:public Delaunay<D,T>{

	PolygonSet<D,T> & polygons;
	std::vector<ConvexPolygon<D,T> > triangulation;

public:
	DelaunayPolypartition(PolygonSet<D,T> & p):polygons(p){

	}

	void update();

	std::vector<ConvexPolygon<D,T> > & getTriangulation(){
		update();
		return triangulation;
	}




private:

};
template<>
DelaunayPolypartition<2,double>::DelaunayPolypartition(PolygonSet<2,double> & p);
template<>
DelaunayPolypartition<2,mpq_class>::DelaunayPolypartition(PolygonSet<2,mpq_class> & p);
template<>
void DelaunayPolypartition<2,double>::update();
template<>
void DelaunayPolypartition<2,mpq_class>::update();
#endif /* DELAUNAYPOLYPARTITION_H_ */
