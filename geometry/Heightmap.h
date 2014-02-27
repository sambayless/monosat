/*
 * ConvexHull.h
 *
 *  Created on: 2014-02-23
 *      Author: sam
 */

#ifndef HEIGHTMAP_H_
#define HEIGHTMAP_H_
#include "GeometryTypes.h"
#include "Polygon.h"


namespace Minisat{

template<unsigned int D,class T=double>
class Heightmap{
public:
	virtual void update()=0;
	virtual Polygon<D,T> & getHeightmesh()=0;
	virtual ~Heightmap(){

	}
};

};
#endif /* CONVEXHULL_H_ */
