/*
 * ConvexHull.h
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#ifndef CONVEXHULL_H_
#define CONVEXHULL_H_
#include "GeometryTypes.h"
#include "ConvexPolygon.h"


namespace Minisat{

template<unsigned int D,class T=double>
class ConvexHull{
public:
	virtual void update()=0;
	virtual ConvexPolygon<D,T> & getHull()=0;
	virtual ~ConvexHull(){

	}
};

};
#endif /* CONVEXHULL_H_ */
