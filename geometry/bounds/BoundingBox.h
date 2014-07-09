/*
 * BoundingBox.h
 *
 *  Created on: Jul 6, 2014
 *      Author: sam
 */

#ifndef BOUNDINGBOX_H_
#define BOUNDINGBOX_H_


#include "BoundingVolume.h"
#include "../Polygon.h"
#include <vector>

//axis-aligned bounding box
template<unsigned int D,class T, class Bound>
class BoundingBox:public BoundingVolume<D,T>{
	Point<D,T> circleCenter;
	T circleRadius;
	Bound & toBound;
public:
	BoundingBox(Bound & toBound):toBound(toBound){

	}
	ShapeType getType(){
		return BOUNDING_BOX;
	}
	void update();
	bool contains(const Point<D,T> & point);
	bool intersects(Shape<D,T> & s);
};
template<unsigned int D,class T>
class BoundingBox<D,T,Polygon<D,T>>:public BoundingVolume<D,T>{
	Point<D,T> max_point;
	Point<D,T> min_point;
	Polygon<D,T> & toBound;
public:
	BoundingBox(Polygon<D,T> & toBound):toBound(toBound){

	}
	ShapeType getType(){
		return BOUNDING_BOX;
	}

	void update(){

		for (int i = 0;i<D;i++){
			max_point[i]=-numeric<T>::infinity();
			min_point[i]=numeric<T>::infinity();
		}
		for (auto & p:toBound){
			for (int i = 0;i<D;i++){
				if(p[i]>max_point[i])
					max_point[i]=p[i];
				if(p[i]<min_point[i])
					min_point[i]=p[i];
			}
		}
	}
	bool contains(const Point<D,T> & point){
		assert(dbg_uptodate());
		for (int i = 0;i<D;i++){
			if(point[i]>max_point[i])
				return false;
			if(point[i]<min_point[i])
				return false;
		}
		return true;
	}
	bool intersects(Shape<D,T> & s){
		assert(false);//not yet implemented
		return false;
	}

private:
	bool dbg_uptodate(){
#ifndef NDEBUG
		//std::vector<Point<D,T>>& vertices  = toBound.getVertices();
		Point<D,T> dbg_max;
		Point<D,T> dbg_min;
		for (int i = 0;i<D;i++){
			dbg_max[i]=-numeric<T>::infinity();
			dbg_min[i]=numeric<T>::infinity();
		}
		for (auto & p:toBound){
			for (int i = 0;i<D;i++){
				if(p[i]>dbg_max[i])
					dbg_max[i]=p[i];
				if(p[i]<dbg_min[i])
					dbg_min[i]=p[i];
			}
		}

		assert(max_point== dbg_max);
		assert(min_point== dbg_min);
#endif
		return true;
	}
};



#endif /* BOUNDINGBOX_H_ */
