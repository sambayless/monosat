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
#include "../ConvexPolygon.h"
#include <vector>

//axis-aligned bounding box

template<unsigned int D,class T>
class AbstractBoundingBox:public BoundingVolume<D,T>{
protected:
	Point<D,T> max_point;
	Point<D,T> min_point;
	//These are important to pre-calculate for the case where we are using expensive arbitrary precision arithmetic.
	Point<D,T> top_left_point;
	Point<D,T> bottom_right_point;

public:
	AbstractBoundingBox(){

	}
	virtual ~AbstractBoundingBox(){

	}
	ShapeType getType(){
		return BOUNDING_BOX;
	}
	virtual void update()=0;
	bool contains(const Point<D,T> & point, bool inclusive=true){

		if(inclusive){
			for (int i = 0;i<D;i++){
				if(point[i]>max_point[i])
					return false;
				if(point[i]<min_point[i])
					return false;
			}
		}else{
			for (int i = 0;i<D;i++){
				if(point[i]>=max_point[i])
					return false;
				if(point[i]<=min_point[i])
					return false;
			}
		}
		return true;
	}
	bool intersects(Shape<D,T> & s, bool inclusive=true){

		if(s.getType()==LINE_SEGMENT){
			if(D==2){
				return intersectsLine2d((LineSegment<2,T>&)s,inclusive );
			}else{
				assert(false);
			}
		}else if (s.getType()==POLYGON || s.getType()== CONVEX_POLYGON ){
			Polygon<D,T> & poly = (Polygon<D,T> &) s;
			if(poly.hasBound()){
				return intersects(*poly.getBound(), inclusive);
			}
			return true;
		}else if (s.getType()== BOUNDING_BOX){
			AbstractBoundingBox<D,T> & B = (AbstractBoundingBox<D,T> & )s;

			for (int i = 0;i<D;i++){
				if(max_point[i]<B.min_point[i] || min_point[i]>B.max_point[i] )
					return false;
			}

			return true;
		}else if (s.getType()== BOUNDING_SPHERE){

			//T sqdist = squared_distance()

		}

		return true;
	}

	private:

		bool intersectsLine2d(LineSegment<2,T> & line, bool inclusive){
			if(line.a == line.b){
				//this is really a point - and it will cause problems below, because the line has no sides.
				return this->contains(line.a,inclusive);
			}

			//check if the corners of the box are all on the same side of the line segment
			int side = (line.whichSide(max_point));
			if(side==0 && inclusive){
				return true;
			}
			int last_side=side;
			side = (line.whichSide(min_point));
			if(side==0){
				if(inclusive)
					return true;
			}else{
				if (last_side != 0 && side != last_side){
					return true;
				}
				last_side=side;
			}
			side = (line.whichSide(top_left_point));
			if(side==0){
				if(inclusive)
					return true;
			}else{
				if (last_side != 0 && side != last_side){
					return true;
				}
				last_side=side;
			}
			side = (line.whichSide(bottom_right_point));
			if(side==0){
				if(inclusive)
					return true;
			}else{
				if (last_side != 0 && side != last_side){
					return true;
				}
			}
			//is this really the only test we need?

			return false;
		}
		bool intersectsConvex2d(ConvexPolygon<2,T> & polygon, bool inclusive){

			return false;
		}
};

template<unsigned int D,class T, class Bound>
class BoundingBox:public AbstractBoundingBox<D,T>{
	Bound & toBound;
public:
	BoundingBox(Bound & toBound): toBound(toBound){

	}
};

template<unsigned int D,class T>
class BoundingBox<D,T,Polygon<D,T>>:public AbstractBoundingBox<D,T>{
	Polygon<D,T> & toBound;
public:
	BoundingBox(Polygon<D,T> & toBound):toBound(toBound){

	}
	ShapeType getType(){
		return BOUNDING_BOX;
	}

	void update(){
		for (int i = 0;i<D;i++){
			this->max_point[i]=-numeric<T>::infinity();
			this->min_point[i]=numeric<T>::infinity();
		}
		for (auto & p:toBound){
			for (int i = 0;i<D;i++){
				if(p[i]>this->max_point[i])
					this->max_point[i]=p[i];
				if(p[i]<this->min_point[i])
					this->min_point[i]=p[i];
			}
		}
		//need to pre-allocate these to make checks using expensive arbitrary precision arithemetic cheap.
		this->top_left_point[0]=this->min_point[0];
		this->top_left_point[1]=this->max_point[1];

		this->bottom_right_point[0]=this->max_point[0];
		this->bottom_right_point[1]=this->min_point[1];
	}

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

		assert(this->max_point== dbg_max);
		assert(this->min_point== dbg_min);
#endif
		return true;
	}
};



#endif /* BOUNDINGBOX_H_ */
