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

#ifndef SHAPE_H_
#define SHAPE_H_
#include <vector>
#include "GeometryTypes.h"

enum ShapeType {
	CONVEX_POLYGON, POLYGON, SHAPE, PLANE, LINE, LINE_SEGMENT, BOUNDING_BOX, BOUNDING_SPHERE
};
/**
 * A concrete shape
 */
template<unsigned int D, class T>
class Shape {
	
public:
	
	Shape() {
	}
	virtual ~Shape() {
	}
	;
	int dimension() {
		return D;
	}
	
	virtual ShapeType getType() {
		return SHAPE;
	}
	
	virtual bool contains(const Point<D, T> & point, bool inclusive)=0;
	virtual bool intersects(Shape<D, T> & s, bool inclusive)=0;
	
};
template<unsigned int D, class T>
std::ostream & operator<<(std::ostream & str, Shape<D, T> & shape) {
	str << "Shape=";
	str << shape.getType();
	str << "";
	return str;
}

#endif /* SHAPE_H_ */
