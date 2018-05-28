/**************************************************************************************************
 The MIT License (MIT)

 Copyright (c) 2014, Sam Bayless

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
 associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
 NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
 OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef BOUNDINGSPHERE_H_
#define BOUNDINGSPHERE_H_
#include "../Polygon.h"
#include <vector>

template<unsigned int D, class T>
class AbstractBoundingSphere: public BoundingVolume<D, T> {
protected:
	Point<D, T> circleCenter;
	T circleRadius;
public:
	AbstractBoundingSphere() {
		
	}
	virtual ~AbstractBoundingSphere() {
		
	}
	ShapeType getType() {
		return BOUNDING_SPHERE;
	}
	virtual void update()=0;
	bool contains(const Point<D, T> & point, bool inclusive = true) {
		if (inclusive) {
			return point.distance_underapprox(circleCenter) <= circleRadius;
		} else {
			return point.distance_underapprox(circleCenter) < circleRadius;
		}
	}
	bool intersects(Shape<D, T> & s, bool inclusive = true) {
		assert(false); //not yet implemented
		return false;
	}
};

template<unsigned int D, class T, class Bound>
class BoundingSphere: public AbstractBoundingSphere<D, T> {
	
	Bound & toBound;
public:
	BoundingSphere(Bound & toBound) :
			toBound(toBound) {
		
	}
	ShapeType getType() {
		return BOUNDING_SPHERE;
	}
	void update();
	
};
template<unsigned int D, class T>
class BoundingSphere<D, T, Polygon<D, T>> : public AbstractBoundingSphere<D, T> {
	Point<D, T> circleCenter;
	T circleRadius;
	Polygon<D, T> & toBound;
public:
	BoundingSphere(Bound & toBound) :
			toBound(toBound) {
		
	}
	ShapeType getType() {
		return BOUNDING_SPHERE;
	}
	void update() {
		circleCenter.zero();
		std::vector<Point<D, T>> & vertices = toBound.getVertices();
		
		for (int i = 0; i < vertices.size(); i++) {
			circleCenter += vertices[i];
		}
		circleCenter /= T(vertices.size());
		circleRadius = T(0);
		for (int i = 0; i < vertices.size(); i++) {
			T dist = circleCenter.distance(vertices[i]);
			if (dist > circleRadius) {
				circleRadius = dist;
			}
		}
	}
	
};

#endif /* BOUNDINGSPHERE_H_ */
