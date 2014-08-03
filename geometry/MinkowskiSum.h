/*
 * MinkowskiSum.h
 *
 *  Created on: Aug 2, 2014
 *      Author: sam
 */

#ifndef MINKOWSKISUM_H_
#define MINKOWSKISUM_H_

#include "GeometryTypes.h"
#include "ConvexPolygon.h"

template<unsigned int D,class T=double>
class MinkowskiSum{
public:
	virtual void update()=0;
	virtual ConvexPolygon<D,T> & getSum()=0;

	virtual ~MinkowskiSum(){

	}
};





#endif /* MINKOWSKISUM_H_ */

