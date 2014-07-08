/*
 * BoundingVolume.h
 *
 *  Created on: Jul 6, 2014
 *      Author: sam
 */

#ifndef BOUNDINGVOLUME_H_
#define BOUNDINGVOLUME_H_

//axis-aligned bounding box
template<unsigned int D,class T>
class BoundingVolume:public Shape<D,T>{
public:
	BoundingVolume(){

	}
	virtual ~BoundingVolume(){}
	virtual void update()=0;
	virtual bool contains(const Point<D,T> & point)=0;
	virtual bool intersects(Shape<D,T> & s)=0;
};



#endif /* BOUNDINGVOLUME_H_ */
