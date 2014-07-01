#ifndef DELAUNAY

#include "GeometryTypes.h"
#include "ConvexPolygon.h"


#include <vector>

template<unsigned int D,class T=double>
class Delaunay{
public:
	virtual void update()=0;
	virtual std::vector<ConvexPolygon<D,T> > & getTriangulation()=0;
	//virtual T getArea()=0;
	virtual ~Delaunay(){

	}
};



#endif
